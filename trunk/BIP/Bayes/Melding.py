# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
# Name:        Melding.py
# Purpose:     The Bayesian melding Class provides
#              uncertainty analyses for simulation models.
#
# Author:      Flávio Codeço Coelho
#
# Created:     2003/08/10
# Copyright:   (c) 2003-2010 by the Author
# Licence:     GPL v3
#-----------------------------------------------------------------------------
import cPickle as CP
import copy
import os
import sqlite3
import sys
import xmlrpclib
from multiprocessing import Pool
from time import time

import numpy
import pylab as P
from numpy import array, nan_to_num, zeros, ones, mean, var, sqrt, floor, isnan,  nansum, median
from numpy.core.records import recarray
from numpy.random import randint, random, seed
from scipy import stats, optimize as optim
from scipy.stats import nanmean
from scipy.stats.kde import gaussian_kde

import PlotMeld as PM
import lhs
import like
from BIP.Bayes.Samplers import MCMC
import pdb


#try:
#    import psyco
#    psyco.full()
#except:
#    pass
try:
    from BIP.Viz.realtime import RTplot
    from liveplots.xmlrpcserver import rpc_plot
    Viz=True
except:
    Viz=False
    print r"""Please install Gnuplot-py to enable realtime visualization.
    http://gnuplot-py.sourceforge.net/
    """

if Viz:
    dtplot = RTplot();phiplot = RTplot();thplot = RTplot()

__docformat__ = "restructuredtext en"

class FitModel(object):
    """
    Fit a model to data generating
    Bayesian posterior distributions of input and
    outputs of the model.
    """
    def __init__(self, K,model, inits,tf,thetanames,phinames,wl=None ,nw=1,verbose=False,burnin=1000, constraints=[]):
        """
        Initialize the model fitter.

        :Parameters:
            - `K`: Number of samples from the priors. On MCMC also the number of samples of the posterior.
            - `model`: Callable (function) returning the output of the model, from a set of parameter values received as argument.
            - `inits`: inits initial values for the model's variables.
            - `tf`: Length of the simulation, in units of time.
            - `phinames`: List of names (strings) with names of the model's variables
            - `thetanames`: List of names (strings) with names of parameters included on the inference.
            - `wl`: window lenght length of the inference window.
            - `nw`: Number of windows to analyze on iterative inference mode
            - `verbose`: Verbosity level: 0, 1 or 2.
            - `burnin`: number of burnin samples, used in the case on mcmc method.
        """
        try:
            assert wl<=tf
        except AssertionError:
           sys.exit("Window Length cannot be larger that Length of the simulation(tf)" )
        assert isinstance(constraints, list)
        self.K = K
        self.L = .1*K if K>2000 else K
        self.finits = inits #first initial values
        self.ftf = tf
        self.full_len =  wl*nw if wl !=None else tf
        self.inits = inits
        self.tf = tf
        self.ew = 0 #expanding windows?
        self.totpop = sum(inits)
        self.model = model
        self.nphi = len(phinames)
        self.ntheta = len(thetanames)
        self.phinames = phinames
        self.thetanames = thetanames
        self.model.func_globals['inits'] = self.inits
        self.model.func_globals['tf'] = self.tf
        self.model.func_globals['thetanames'] = self.thetanames
        self.wl = wl
        self.nw = nw
        self.done_running = False
        self.prior_set = False
        self.burnin = burnin
        self.verbose = verbose
        self.constraints = constraints
        self.pool = False #this will be updated by the run method.
        self.Me = Meld(K=self.K,L=self.L,model=self.model,ntheta=self.ntheta,nphi=self.nphi,verbose=self.verbose)
        self.AIC = 0
        self.BIC = 0
        self.DIC = 0
        # To be defined by self.set.priors
        self.pdists = None
        self.ppars = None
        self.plims = None
        self.tdists = None
        self.tpars = None
        self.tlims = None

    def _plot_MAP(self,data,pmap):
        """
        Generates a plot of a full run of the model parameterized with the maximum a posteriori
        estimates of the parameters.

        :Parameters:
            - `data`: data dictionary as passed to optimize
            - `pmap`: MAP parameter values
        """
        p = RTplot(persist=1)
        self.model.func_globals['inits'] = self.finits; self.model.func_globals['tf'] = self.full_len
        simseries = self.model(list(pmap))
        self.model.func_globals['inits'] = self.inits; self.model.func_globals['tf'] = self.tf
        if 'time' in data:
            data.pop('time')
        for n,d in data.items():
            p.plotlines(nan_to_num(d).tolist(),range(len(d)),['Obs. %s'%n], '','points', 0)
            v=self.phinames.index(n)
            p.plotlines(simseries.T[v].tolist(),range(len(d)), data.keys(), "Simulation with MAP parameters %s=%s"%(self.thetanames,pmap))

    def AIC_from_RSS(self,):
        """
        Calculates the Akaike information criterion from the residual sum of squares 
        of the best fitting run.
        """
        pass

    def optimize(self, data, p0, optimizer='scipy', tol=0.0001, verbose=0, plot=0):
        """
        Finds best parameters using an optimization approach

        :Parameters:
            - `data`: Dictionary of observed series
            - `p0`: Sequence (list or tuple) of initial values for the parameters
            - `optimizer`: Optimization library to use: 'scipy': fmin (Nelder-Mead) or 'oo':OpenOpt.NLP
            - `tol`: Tolerance of the error
            - `verbose`: If true show stats of the optimization run at the end
            - `plot`: If true plots a run based on the optimized parameters.
        """
        try:
            import openopt
            oo = True
        except ImportError:
            oo = False
        assert isinstance(p0,(list,tuple))
        def mse(theta):
            s1 = model_as_ra(theta,  self.model,  self.phinames)
            return self._rms_error(s1, data)
        if optimizer == "scipy":
            potimo = optim.fmin(mse,p0,ftol=tol, disp=verbose)
        elif optimizer == "oo":
            if not oo:
                potimo = optim.fmin(mse,p0,ftol=tol, disp=verbose)
            else: 
                lb = array([l[0] for l in self.tlims]) if self.tlims else array([-numpy.inf]*self.ntheta)
                ub = array([l[1] for l in self.tlims]) if self.tlims else array([numpy.inf]*self.ntheta)
                p = openopt.NLP(mse, p0, lb=lb, ub=ub, ftol=tol, iprint=10)
                p.solve('ralg')
                potimo = p.xf
        else:#use fmin as fallback method
            potimo = optim.fmin(mse,p0,ftol=tol, disp=verbose)
        if plot:
            self._plot_MAP(data, potimo)
        return potimo
        
    def _rms_error(self, s1, s2):
        '''
        Calculates a the error between a model-
        generated time series and a observed time series.
        It uses a normalized RMS deviation.

        :Parameters:
            - `s1`: model-generated time series. 
            - `s2`: observed time series. dictionary with keys matching names of s1
        :Types:
            - `s1`: Record array or list.
            - `s2`: Dictionary or list; must be a dictionary if s1 is a RA
        
        s1 and s2 can also be both lists of lists or lists of arrays of the same length.

        :Return:
            The Root mean square deviation between `s1` and `s2`.
        '''
        if isinstance(s1, recarray):
            assert isinstance(s2, dict)
            err = []
            for k in s2.keys():
                if k not in s1.dtype.names:
                    continue
                ls1 = len(s1[k]) #handles the cases where data is slightly longer that simulated series.
#                print s2[k]
#                try:
#                pdb.set_trace()
                if len(s2[k].shape) >1:
                    s2[k] = s2[k].mean(axis=1)
                dif = s1[k]-s2[k][:ls1].astype(float)
                
                dif[isnan(dif)] = 0
                e = sqrt(mean(dif**2))
#                except TypeError:
#                    print s1[k], s2[k]

                err.append(e) 
        elif isinstance(s1, list):
            assert isinstance(s2, list) and len(s1) ==len(s2)
            s1 = array(s1).astype(float)
            s2 = array(s2).astype(float)
            err = [sqrt(nanmean((s-t)**2)) for s, t in zip(s1, s2) ]
            #err = [sum((s-t)**2./t**2) for s, t in zip(s1, s2)]
        rmsd = nan_to_num(mean(err))
        return rmsd
        
    def set_priors(self,tdists,tpars, tlims,pdists,ppars,plims):
        """
        Set the prior distributions for Phi and Theta

        :Parameters:
            - `pdists`: distributions for the output variables. For example: [scipy.stats.uniform,scipy.stats.norm]
            - `ppars`: paramenters for the distributions in pdists. For example: [(0,1),(0,1)]
            - `plims`: Limits of the range of each phi. List of (min,max) tuples.
            - `tdists`: same as pdists, but for input parameters (Theta).
            - `tpars`: same as ppars, but for tdists.
            - `tlims`: Limits of the range of each theta. List of (min,max) tuples.
        """
        self.pdists = pdists
        self.ppars = ppars
        self.plims = plims
        self.tdists = tdists
        self.tpars = tpars
        self.tlims = tlims
        self._init_priors()
        self.prior_set = True
    
    def prior_sample(self):
        """
        Generates a set of sample from the starting theta prior distributions
        for reporting purposes.
        
        :Returns:
            Dictionary with (name,sample) pairs
        """
        s = {}
        for i, n in enumerate(self.thetanames):
            s[n]=self.tdists[i](*self.tpars[i]).rvs(self.K)
        return s
            
    def _init_priors(self, prior=None):
        """
        Initialize priors either from distributions or previous posteriors
        """
        if prior!=None and prior['theta'] and prior['phi']:
            self.Me.setThetaFromData(self.thetanames,prior['theta'],self.tlims)
            self.Me.setPhiFromData(self.phinames,prior['phi'],self.plims)
        else:
            print "++++>"
            self.Me.setTheta(self.thetanames,self.tdists,self.tpars, self.tlims)
            self.Me.setPhi(self.phinames,self.pdists,self.ppars,self.plims)
            
    def do_inference(self, prior, data, predlen, method, likvar):
        """
        """
        self._init_priors(prior)
        succ=0
        att = 1
        for n in data.keys():
            if n not in self.phinames:
                data.pop(n)
        if method == "SIR":
            while not succ: #run sir Until is able to get a fit
                print 'attempt #',att
                succ = self.Me.sir(data=data,variance=likvar,pool=self.pool,t=self.tf)
                att += 1
            pt,series = self.Me.getPosteriors(t=self.tf)
        elif method == "MCMC":
            while not succ: #run sir Until is able to get a fitd == "mcmc":
                print 'attempt #',att
                succ = self.Me.mcmc_run(data,t=self.tf,likvariance=likvar,burnin=self.burnin, method = 'MH')
            pt = self.Me.post_theta
            series = self.Me.post_phi
        elif method == "DREAM":
            while not succ: #run sir Until is able to get a fitd == "mcmc":
                print 'attempt #',att
                succ = self.Me.mcmc_run(data,t=self.tf,likvariance=likvar,burnin=self.burnin,  method='dream', constraints=self.constraints)
            pt = self.Me.post_theta
            series = self.Me.post_phi
        elif method == "ABC":
            #TODO: allow passing of fitfun
            self.Me.abcRun(data=data,fitfun=None,pool=self.pool, t=self.tf)
            pt,series = self.Me.getPosteriors(t=self.tf)
        if self.Me.stop_now:
            # if fitting has been prematurely interrupted by user
            return None,None,None, None, None
        pp = series[:,-1]
        # TODO: figure out what to do by default with inits
        if self.nw >1 and self.adjinits and not self.ew:
            adiff = array([abs(pp[vn]-data[vn][-1]) for vn in data.keys()])
            diff = adiff.sum(axis=0) #sum errors for all oserved variables
            initind = diff.tolist().index(min(diff))
            self.inits = [pp[vn][initind] for vn in self.phinames]
            for i, v in enumerate(self.phinames):
                if v in data.keys():
                    self.inits[i] = data[v][-1] 
            self.model.func_globals['inits'] = self.inits
        
        if predlen:
            predseries = self.Me.getPosteriors(predlen)[1]
        return pt,pp, series,predseries,att

    def _save_to_db(self, dbname, data):
        '''
        Saves data to a sqlite3 db.
        
        :Parameters:
            - `dbname`: name of the database file
            - `data`: Data dictionary as created by `format_db_tables` method.
        '''
        def row_generator(var):
            '''
            var is a dictionary.
            '''
#            if isinstance(var, numpy.recarray):
#                for repl in var:
#                    if isinstance(repl, numpy.recarray): # this is the case of variables, where each replicate is a time series
#                        for t in repl:
#                            if not isinstance(t, numpy.core.records.record):
#                                t=(t,)
#                            yield tuple(t) #variable tuple at time t in replicate repl
#                    else: #this is the case of parameters
#                        yield tuple(repl)
#            elif isinstance(var, dict):
#                try:
            for r in zip(*var.values()):
                if not isinstance(r, tuple):
                    r = (r,)
                t = []
                for i in r:
                    if isinstance(i, numpy.ndarray):
                        t+=i.tolist()
                    else:
                        t += [i]
                    r = tuple(t)
                yield r
#                except:
#                    print var.keys(),  var.values()
                    
        create = True
        if not dbname.endswith('.sqlite'):
            dbname += '.sqlite'
        if os.path.exists(dbname):
            create = False
        con = sqlite3.connect(dbname)
        #create tables
        for k, v in data.items():
            if isinstance(v, dict):
                labs = []
                for k2, v2 in v.items():
                    if isinstance(v2, numpy.ndarray) and len(v2.shape)>1:
                        labs+=[k2+str(i) for i in range(v2.shape[1])]
                    else:
                        labs.append(k2)
                nv = len(labs)#+1 #variables plus primary key
                tstrc = k+'(pk integer primary key asc autoincrement,'+','.join(labs)+')'
                tstr = k+'('+','.join(labs)+')'
                if create:
                    con.execute('create table '+tstrc)
#            elif isinstance(v, numpy.recarray):
#                nv = len(v.dtype.names) +1 #variables plus primary key
#                tstr = k+"("+','.join(v.dtype.names)+')'
#                if create:
#                    con.execute("create table "+ tstrc)
            else:
                raise TypeError("Non-valid data structure.")
            #print "insert into "+tstr+" values("+",".join(['?']*nv)+")"
            con.executemany("insert into "+tstr+" values("+",".join(['?']*nv)+")", row_generator(v))
        con.commit()
        con.close()
    def run(self, data,method,likvar,pool=False,adjinits=True,ew=0, dbname='results', monitor=False, initheta=[]):
        """
        Fit the model against data

        :Parameters:
            - `data`: dictionary with variable names and observed series, as Key and value respectively.
            - `method`: Inference method: "ABC", "SIR", "MCMC" or "DREAM"
            - `likvar`: Variance of the likelihood function in the SIR and MCMC method
            - `pool`: Pool priors on model's outputs.
            - `adjinits`: whether to adjust inits to data
            - `ew`: Whether to use expanding windows instead of moving ones.
            - `dbname`: name of the sqlite3 database
            - `monitor`: Whether to monitor certains variables during the inference. If not False, should be a list of valid phi variable names.
            - `initheta`: starting position in parameter space for the sampling to start. (only used by MCMC and DREAM) 
        """
        self.ew = ew
        self.adjinits = adjinits
        self.pool = pool
        assert isinstance(initheta, list) #type checking
        self.Me.initheta = initheta
        if not self.prior_set: return
        if monitor:
            self._monitor_setup()
        start = time()
        d = data
        prior = {'theta':[],'phi':[]}
        os.system('rm %s_*.pickle'%dbname)
        if self.wl == None:
            self.wl = floor(len(d.values()[0])/self.nw)
        wl = self.wl
        for w in range(self.nw):
            t0 = time()
            print '==> Window # %s of %s!'%(w+1,self.nw)
            if w>0:
                rt = (self.nw-w+1)*tel
                print "==> Remaining time: %s minutes and %s seconds."%(rt//60, rt%60)
            self.tf=(w+1)*wl if ew else wl #setting tf according to window type
            self.model.func_globals['tf'] = self.tf
            d2 = {}
            for k,v in d.items():#Slicing data to the current window
                d2[k] = v[:self.tf] if ew else v[w*wl:w*wl+wl]            
            if w==0 and adjinits:
                for n in d2.keys():
                    if n not in self.phinames:
                        continue
                    i = self.phinames.index(n)
                    self.inits[i] = nan_to_num(d2[n][0]) if  not isinstance(d2[n][0], numpy.ndarray ) else nan_to_num(nanmean(d2[n][0]))
                    #TODO: figure out how to balance the total pop
#                    self.inits[0] += self.totpop-sum(self.inits) #adjusting susceptibles
                    self.model.func_globals['inits'] = self.inits
            pt,pp, series,predseries,att = self.do_inference(data=d2, prior=prior,predlen=wl, method=method,likvar=likvar)
            if self.Me.stop_now:
                return
            self.AIC += 2. * (self.ntheta - self.Me.likmax) # 2k - 2 ln(L)
            self.BIC += self.ntheta * numpy.log(self.wl*len(d2)) - 2. * self.Me.likmax # k ln(n) - 2 ln(L)
            self.DIC = self.Me.DIC
            # ===Saving results===
            with open('%s_%s%s'%(dbname, w, ".pickle"),'w') as f:
                #save weekly posteriors of theta and phi, posteriors of series, data (d) and predictions(z)
                CP.dump((pt,series,d,predseries, att*self.K),f)
            if dbname:
                if os.path.exists(dbname+".sqlite") and w ==0:
                    os.remove(dbname+".sqlite")
                self._format_db_tables(dbname, w, data, pt, series, predseries, self.AIC, self.BIC, self.DIC)
            prior = {'theta':[],'phi':[]}
            for n in pt.dtype.names:
                prior['theta'].append(pt[n])
            #beta,alpha,sigma,Ri  = median(pt.beta),median(pt.alpha),median(pt.sigma),median(pt.Ri
            for n in pp.dtype.names:
                #print compress(isinf(pp[n]),pp[n])
                prior['phi'].append(pp[n])
            if monitor:
                self._monitor_plot(series,prior,d2,w,data,vars=monitor)
            
            tel = time()-t0
        self.AIC /=self.nw #averaging
        self.BIC /=self.nw
        self.Me.AIC = self.AIC
        self.Me.BIC = self.BIC
        self.DIC = self.Me.DIC
        print "time: %s seconds"%(time()-start)
        
        self.done_running = True
    
    def _format_db_tables(self, dbname, w, data, pt, series, predseries, AIC, BIC, DIC):
        """
        Formats results for writing to database
        """
        #TODO: Write tests for this
        # data table does not require formatting
        # Dates for this window
        self.wl=int(self.wl)
        if 'time' in data:
            if isinstance(data['time'],  numpy.ndarray):
                ts = data['time'].tolist()
            else:
                ts = data['time']
            dates = ts[w*self.wl:w*self.wl+self.wl]
            preddates = ts[(w+1)*self.wl:(w+1)*self.wl+self.wl]
        else:
            dates = range(w*self.wl, w*self.wl+self.wl)
            preddates = range((w+1)*self.wl, (w+1)*self.wl+self.wl)
        # Parameters table
        ptd = {'time':[dates[-1]]*len(pt[pt.dtype.names[0]])}
        for n in pt.dtype.names:
            ptd[n] = pt[n]
        # Series table
        seriesd = {'time':dates*series.shape[0]}
        predseriesd = {'time':preddates*series.shape[0]}
        for n in series.dtype.names:
            seriesd[n] = series[n].ravel()
        for n in predseries.dtype.names:
            predseriesd[n] = predseries[n].ravel()
        # AIC and BIC table
        gof = {'time':[dates[-1]], 'AIC':[AIC], 'BIC':[BIC], 'DIC':[DIC]}
        self._save_to_db(dbname, {'post_theta':ptd, 
                                                                     'series':seriesd, 
                                                                     'data': data, 
                                                                     'predicted_series':predseriesd, 
                                                                     'GOF':gof,
                                                                     })
    
    def _monitor_setup(self):
        """
        Sets up realtime plotting for inference
        """
        #theta histograms (for current window)
        self.hst = xmlrpclib.ServerProxy('http://localhost:%s'%rpc_plot(hold=1), allow_none=1)#RTplot() 
        #full data and simulated series
        self.fsp = xmlrpclib.ServerProxy('http://localhost:%s'%rpc_plot(hold=1), allow_none=1)#RTplot()
        # phi time series (model output for the current window)
        self.ser = xmlrpclib.ServerProxy('http://localhost:%s'%rpc_plot(hold=1), allow_none=1)#RTplot()

    def _get95_bands(self,series,vname):
        '''
        Returns 95% bands for series of all variables in vname
        
        :Parameters:
            - `series`: record array containing the series
            - `vname`: list of strings of variables in series for which we want to get the bands
        '''
        i5 = array([stats.scoreatpercentile(t,2.5) for t in series[vname].T])
        i95 = array([stats.scoreatpercentile(t,97.5) for t in series[vname].T])
        return i5,i95
        
    def _long_term_prediction_plot(self, cpars,vind,  w):
        """
        Plots the simulated trajectory predicted from best fit parameters.
        
        :Parameters:
            - `cpars`: best fit  parameter set
            - `vind`: List with indices(in self.phinames) to variables to be plotted 
            - `w`: current window number.
        """
        if self.full_len-(self.wl*(w+1)) == 0:
            return
        self.model.func_globals['tf'] = self.full_len if self.ew else self.full_len-(self.wl*(w+1))
        if self.ew: self.model.func_globals['inits'] = self.finits
        simseries = self.model(cpars)
        simseries = [simseries[:, i].tolist() for i in range(self.nphi) if i in vind]
        snames = [n for n in self.phinames if i in vind]
        self.model.func_globals['tf'] = self.tf
        xinit = 0 if self.ew else self.wl*w+self.wl
#        print xinit, xinit+len(simseries[0])
        self.fsp.lines(simseries,range(xinit,xinit+len(simseries[0])), snames, "Best fit simulation after window %s"%(w+1))
    
        
    
    def _monitor_plot(self, series, prior, d2,w,data, vars):
        """
        Plots real time data
        
        :Parameters:
            - `series`: Record array with the simulated series.
            - `prior`: Dctionary with the prior sample of Theta
            - `d2`: Dictionary with data for the current fitting window.
            - `w`: Integer id of the current fitting window.
            - `data`: Dictionary with the full dataset.
            - `vars`: List with variable names to be plotted.
        """
#        diff = 0
        for vn in d2.keys():
            #sum errors for all oserved variables
            if isinstance(d2[vn], numpy.ndarray) and len(d2[vn].shape)>1:
                diff = abs(series[vn][:, -1]-nanmean(d2[vn][-1]))
            else:
                diff = abs(series[vn][:, -1]-d2[vn][-1])
#        print diff, min(diff)
        initind = diff.tolist().index(min(diff))
        vindices = [self.phinames.index(n) for n in vars]
        for n in vars:
            if n not in d2:
                continue
            i5,i95 = self._get95_bands(series,n)
            self.ser.lines(series[n][initind].tolist(), None,  ["Best run's %s"%n], 'Window %s'%(w+1))
            self.ser.lines(i5.tolist(),None, ['2.5%'], 'Window %s'%(w+1))
            self.ser.lines(i95.tolist(),None,  ['97.5%'],'Window %s'%(w+1))
            self.ser.lines(d2[n].tolist(),None, ['Obs. %s'%n], 'Window %s'%(w+1), 'points')
            self.fsp.lines(data[n].T.tolist(),None, ['Obs. %s'%n], 'Window %s'%(w+1), 'points')
        self.hst.histogram(array(prior['theta']).tolist(), self.thetanames,'Window %s'%(w+1),  1)
        cpars = [prior['theta'][i][initind] for i in range(self.ntheta)]
        self._long_term_prediction_plot(cpars,vindices, w)
        self.ser.clearFig()
        self.hst.clearFig()
        self.fsp.clearFig()
    
    def plot_results(self, names=[],  dbname="results", savefigs=0):
        """
        Plot the final results of the inference
        """
        if not names:
            names = self.phinames
        try: #read the data files
            pt,series,predseries,obs = self._read_results(dbname)
        except:
            if not self.done_running:
                return
        if obs.has_key('time'):
            tim = numpy.array(obs['time'])
        else:
            tim = numpy.arange(self.nw*self.wl)
        #PM.plot_par_series(range(len(pt)),pt)
        priors = self.prior_sample()
        PM.plot_par_violin(tim[self.wl-1::self.wl],pt, priors)
        P.xlabel('windows')
        if savefigs:
            P.savefig(dbname+"_violin.svg")
        P.figure()
        PM.plot_series2(tim,obs,series,names=names, wl=self.wl)
        if savefigs:
                P.savefig(dbname+"_series.svg")
        if self.nw > 1:
            P.figure()
            PM.pred_new_cases(obs,predseries,self.nw,names,self.wl)
            P.gca().legend(loc=0)
            P.xlabel('windows')
            P.figure()
            PM.plot_series2(tim,obs,predseries,names=names,
                            title='Predicted vs. Observed series',lag=True)
            P.xlabel('windows')
            if savefigs:
                P.savefig(dbname+"_predseries.svg")
        P.show()

    def _read_results(self, nam):
        """
        read results from disk
        """
        pt,series,predseries = [],[],[]

        for w in range(self.nw):
            fn = "%s_%s.pickle"%(nam, w)
            print fn
            f = open(fn,'r')
            a,b,obs,pred, samples = CP.load(f)
            f.close()
            pt.append(a)
            series.append(b)
            predseries.append(pred)
        return pt,series,predseries,obs
            

class Meld(object):
    """
    Bayesian Melding class
    """
    def __init__(self,  K,  L, model, ntheta, nphi, alpha = 0.5, verbose = 0, viz=False ):
        """
        Initializes the Melding class.
        
        :Parameters:
            - `K`: Number of replicates of the model run. Also determines the prior sample size.
            - `L`: Number of samples from the Posterior distributions. Usually 10% of K.
            - `model`: Callable taking theta as argument and returning phi = M(theta).
            - `ntheta`: Number of inputs to the model (parameters).
            - `nphi`: Number of outputs of the model (State-variables)
            - `verbose`: 0,1, 2: whether to show more information about the computations
            - `viz`: Boolean. Wether to show graphical outputs of the fitting process
        """
        self.K = K
        self.L = L
        self.verbose = verbose
        self.model = model
        self.likelist = [] #list of likelihoods
        self.q1theta = recarray(K,formats=['f8']*ntheta) #Theta Priors (record array)
        self.post_theta = recarray(L,formats=['f8']*ntheta) #Theta Posteriors (record array)
        self.q2phi = recarray(K,formats=['f8']*nphi) #Phi Priors (record array)
        self.phi = recarray(K,formats=['f8']*nphi) #Phi model-induced Priors (record array)
        self.q2type = [] #list of distribution types
        self.post_phi = recarray(L,formats=['f8']*nphi) #Phi Posteriors (record array)
        self.ntheta = ntheta
        self.nphi = nphi
        self.thetapars = []
        self.phipars = []
        self.alpha = alpha #pooling weight of user-provided phi priors
        self.done_running = False
        self.theta_dists = {}#parameterized rngs for each parameter
        self.phi_dists = {}#parameterized rngs for each variable
        self.likmax = -numpy.inf
        self.initheta = []
        self.AIC = None
        self.BIC = None
        self.DIC = None
        self.stop_now = False #Flag to interrupt fitting loop. inner loop checks this flag and exit when it is true.
        self.current_step = 0
        self.proposal_variance = 0.0000001
        self.adaptscalefactor = 1 #adaptive variance. Used my metropolis hastings
        self.salt_band = 0.1 #relaxation factor of the prior limits
        if Viz: #Gnuplot installed
            self.viz = viz
        else:
            self.viz = False
        if self.verbose == 2:
            self.every_run_plot = xmlrpclib.ServerProxy('http://localhost:%s'%rpc_plot(hold=1), allow_none=1)
        self.po = Pool() #pool of processes for parallel processing
    
    def current_plot(self, series, data, idx, vars=[], step=0):
        """
        Plots the last simulated series
        
         :Parameters:
            - `series`: Record array with the simulated series.
            - `idx`: Integer index of the curve to plot .
            - `data`: Dictionary with the full dataset.
            - `vars`: List with variable names to be plotted.
            - `step`: Step of the chain
        """
        try:
            if self.lastidx == idx:
                return
        except AttributeError:
            pass
        #print series.shape, idx
        best_series = [series[k][idx].tolist() for k in data.keys()]
        d = [data[k].tolist() for k in data.keys()]
        self.every_run_plot.lines(d,[],data.keys(), "Best fit. Last updated on %s"%step, 'points' )
        self.every_run_plot.lines(best_series,[],data.keys(), "Best fit. Last updated on %s"%step, 'lines' )
        self.every_run_plot.clearFig()
        self.lastidx = idx
    
    def setPhi(self, names, dists=[stats.norm], pars=[(0, 1)], limits=[(-5,5)]):
        """
        Setup the models Outputs, or Phi, and generate the samples from prior distributions 
        needed for the melding replicates.
        
        :Parameters:
            - `names`: list of string with the names of the variables.
            - `dists`: is a list of RNG from scipy.stats
            - `pars`: is a list of tuples of variables for each prior distribution, respectively.
            - `limits`: lower and upper limits on the support of variables.
        """
        if len(names) != self.nphi:
            raise ValueError("Number of names(%s) does not match the number of output variables(%s)."%(len(names),self.nphi))
        self.q2phi.dtype.names = names
        self.phi.dtype.names = names
        self.post_phi.dtype.names = names
        self.plimits = limits
        self.phipars = pars
        for n,d,p in zip(names,dists,pars):
            self.q2phi[n] = lhs.lhs(d,p,self.K).ravel()
            self.q2type.append(d.name)
            self.phi_dists[n]=d(*p)


        
    def setTheta(self, names, dists=[stats.norm], pars=[(0, 1)], lims=[(0, 1)]):
        """
        Setup the models inputs and generate the samples from prior distributions 
        needed for the dists the melding replicates.
        
        :Parameters:
            - `names`: list of string with the names of the parameters.
            - `dists`: is a list of RNG from scipy.stats
            - `pars`: is a list of tuples of parameters for each prior distribution, respectivelydists
        """
        self.q1theta.dtype.names = names
        self.post_theta.dtype.names = names
        self.thetapars = pars
        self.tlimits = lims
        if os.path.exists('q1theta'):
            self.q1theta = CP.load(open('q1theta','r'))
        else:
            for n,d,p in zip(names,dists,pars):
                self.q1theta[n] = lhs.lhs(d,p,self.K).ravel()
                self.theta_dists[n]=d(*p)
        
    def add_salt(self,dataset,band):
        """
        Adds a few extra uniformly distributed data 
        points beyond the dataset range.
        This is done by adding from a uniform dist.
        
        :Parameters:
            - `dataset`: vector of data
            - `band`: Fraction of range to extend [0,1[
            
        :Returns:
            Salted dataset.
        """
        dmax = max(dataset)
        dmin = min(dataset)
        drange = dmax-dmin
        hb = drange*band/2.
        d = numpy.concatenate((dataset,stats.uniform(dmin-hb,dmax-dmin+hb).rvs(self.K*.05)))
        return d
    
    def setThetaFromData(self,names,data, limits):
        """
        Setup the model inputs and set the prior distributions from the vectors
        in data.
        This method is to be used when the prior distributions are available in 
        the form of a sample from an empirical distribution such as a bayesian
        posterior.
        In order to expand the samples provided, K samples are generated from a
        kernel density estimate of the original sample.
        
        :Parameters:
            - `names`: list of string with the names of the parameters.
            - `data`: list of vectors. Samples of a proposed distribution
            - `limits`: List of (min,max) tuples for each theta to make sure samples are not generated outside these limits.
        """
        self.q1theta.dtype.names = names
        self.post_theta.dtype.names = names
        self.tlimits = limits
        tlimits = dict(zip(names,limits))
        class Proposal:
            "class wrapping a kde adding similar interface to stats dists"
            def __init__(self,name, dist):
                self.name = name
                self.dist = dist
            def __call__(self):
                smp = -numpy.inf
                while not (smp>=tlimits[self.name][0] and smp<=tlimits[self.name][1]):
                    smp = self.dist.resample(1)[0][0]
                return smp
            def rvs(self, *args, **kwds):
                if 'size' in kwds:
                    sz = kwds['size']
                else:
                    sz = 1
                    #TODO implement return of multiple samples here if necessary
                smp = -numpy.inf
                while not (smp>=tlimits[self.name][0] and smp<=tlimits[self.name][1]):
                    smp = self.dist.resample(sz)[0][0]
                return smp
            def stats(self, moments):
                if isinstance(self.dist, (stats.rv_continuous, stats.rv_discrete)):
                    return self.dist.stats(moments='m')
                else:
                    return self.dist.dataset.mean()
            def pdf(self, x):
                return self.dist.evaluate(x)[0]
            def pmf(self, x):
                return self.dist.evaluate(x)[0]
                
        if os.path.exists('q1theta'):
            self.q1theta = CP.load(open('q1theta','r'))
        else:
            for n,d,lim in zip(names,data,limits):
                smp = []
                #add some points uniformly across the support 
                #to help MCMC to escape from prior bounds
                salted = self.add_salt(d,self.salt_band)

                dist = gaussian_kde(salted)
                while len(smp)<self.K:
                    smp += [x for x in dist.resample(self.K)[0] if x >= lim[0] and x <= lim[1]]
                #print self.q1theta[n].shape, array(smp[:self.K]).shape
                self.q1theta[n] = array(smp[:self.K])
                self.theta_dists[n] = Proposal(copy.deepcopy(n),copy.deepcopy(dist))

    def setPhiFromData(self,names,data,limits):
        """
        Setup the model outputs and set their prior distributions from the
        vectors in data.
        This method is to be used when the prior distributions are available in
        the form of a sample from an empirical distribution such as a bayesian
        posterior.
        In order to expand the samples provided, K samples are generated from a
        kernel density estimate of the original sample.

        :Parameters:
            - `names`: list of string with the names of the variables.
            - `data`: list of vectors. Samples of the proposed distribution.
            - `limits`: list of tuples (ll,ul),lower and upper limits on the support of variables.
        """
        self.q2phi.dtype.names = names
        self.phi.dtype.names = names
        self.post_phi.dtype.names = names
        self.limits = limits
        for n,d in zip(names,data):
            i = 0
            smp = []
            while len(smp)<self.K:
                try:
                    smp += [x for x in gaussian_kde(d).resample(self.K)[0] if x >= limits[i][0] and x <= limits[i][1]]
                except:
                    #d is has no variation, i.e., all elements are the same
                    #print d
                    #raise LinAlgError, "Singular matrix"
                    smp = ones(self.K)*d[0] #in this case return a constant sample
            self.q2phi[n] = array(smp[:self.K])
            self.q2type.append('empirical')
            i += 1
        #self.q2phi = self.filtM(self.q2phi, self.q2phi, limits)

        
    def run(self,*args):
        """
        Runs the model through the Melding inference.model
        model is a callable which return the output of the deterministic model,
        i.e. the model itself.
        The model is run self.K times to obtain phi = M(theta).
        """
        #TODO: implement calculation of AIC,BIC and DIC for SIR
        for i in xrange(self.K):
            theta = [self.q1theta[n][i] for n in self.q1theta.dtype.names]
            r = self.po.apply_async(self.model, theta)
            res = r.get()
            self.phi[i]= res[-1]#self.model(*theta)[-1] #phi is the last point in the simulation

        self.done_running = True
        
    def getPosteriors(self,t):
        """
        Updates the posteriors of the model's output for the last t time steps.
        Returns two record arrays:
        - The posteriors of the Theta
        - the posterior of Phi last t values of time-series. self.L by `t` arrays.

        :Parameters:
            - `t`: length of the posterior time-series to return.
        """
        if not self.done_running:
            print "Estimation has not yet been run"
            return
        if t > 1:
            self.post_phi = recarray((self.L,t),formats=['f8']*self.nphi)
            self.post_phi.dtype.names = self.phi.dtype.names
        def cb(r):
            '''
            callback function for the asynchronous model runs.
            r: tuple with results of simulation (results, run#)
            '''
            if t == 1:
                self.post_phi[r[1]] = (r[0][-1],)
                #self.post_phi[r[1]]= [tuple(l) for l in r[0][-t:]]
            else:
                self.post_phi[r[1]]= [tuple(l) for l in r[0][-t:]]
        po = Pool()
        #random indices for the marginal posteriors of theta
        pti = lhs.lhs(stats.randint,(0,self.L),siz=(self.ntheta,self.L))
        for i in xrange(self.L):#Monte Carlo with values of the posterior of Theta
            theta = [self.post_theta[n][pti[j,i]] for j,n in enumerate(self.post_theta.dtype.names)]
            po.apply_async(enumRun, (self.model,theta,i), callback=cb)
            if i%100 == 0 and self.verbose:
                print "==> L = %s\r"%i

        po.close()
        po.join()
        return self.post_theta, self.post_phi

    def filtM(self,cond,x,limits):
        '''
        Multiple condition filtering.
        Remove values in x[i], if corresponding values in
        cond[i] are less than limits[i][0] or greater than
        limits[i][1].

        :Parameters:
            - `cond`: is an array of conditions.
            - `limits`: is a list of tuples (ll,ul) with length equal to number of lines in `cond` and `x`.
            - `x`: array to be filtered.
        '''
        # Deconstruct the record array, if necessary.
        names = []
        if isinstance(cond, recarray):
            names = list(cond.dtype.names)
            cond = [cond[v] for v in cond.dtype.names]
            x = [x[v] for v in x.dtype.names]

        cond = array(cond)
        cnd = ones(cond.shape[1],int)
        for i,j in zip(cond,limits):
            ll = j[0]
            ul = j[1]
            #print cond.shape,cnd.shape,i.shape,ll,ul
            cnd = cnd & less(i,ul) & greater(i,ll)
        f = compress(cnd,x, axis=1)

        if names:#Reconstruct the record array
            r = recarray((1,f.shape[1]),formats=['f8']*len(names),names=names)
            for i,n in enumerate(names):
                r[n]=f[i]
            f=r

        return f

        

    def logPooling(self,phi):
        """
        Returns the probability associated with each phi[i]
        on the pooled pdf of phi and q2phi.

        :Parameters:
            - `phi`: prior of Phi induced by the model and q1theta.
        """
        
#       Estimating the multivariate joint probability densities
        phidens = gaussian_kde(array([phi[n][:,-1] for n in phi.dtype.names]))

        q2dens = gaussian_kde(array([self.q2phi[n] for n in self.q2phi.dtype.names]))
#       Determining the pooled probabilities for each phi[i]
#        qtilphi = zeros(self.K)
        lastp = array([list(phi[i,-1]) for i in xrange(self.K)])
#        print lastp,lastp.shape
        qtilphi = (phidens.evaluate(lastp.T)**(1-self.alpha))*q2dens.evaluate(lastp.T)**self.alpha
        return qtilphi/sum(qtilphi)

    def abcRun(self,fitfun=None, data={}, t=1,pool=False,savetemp=False):
        """
        Runs the model for inference through Approximate Bayes Computation
        techniques. This method should be used as an alternative to the sir.
        
        :Parameters:
             - `fitfun`: Callable which will return the goodness of fit of the model to data as a number between 0-1, with 1 meaning perfect fit
             - `t`: number of time steps to retain at the end of the of the model run for fitting purposes.
             - `data`: dict containing observed time series (lists of length t) of the state variables. This dict must have as many items the number of state variables, with labels matching variables names. Unorbserved variables must have an empty list as value.
             - `pool`: if True, Pools the user provided priors on the model's outputs, with the model induced priors.
             - `savetemp`: Should temp results be saved. Useful for long runs. Alows for resuming the simulation from last sa
        """
        seed()
        if not fitfun:
            fitfun = basicfit
        if savetemp:
            CP.dump(self.q1theta,open('q1theta','w'))
#       Running the model ==========================
        phi = self.runModel(savetemp,t)

        print "==> Done Running the K replicates\n"
        # Do Log Pooling
        if not pool:
            qtilphi = ones(self.K)
        else:
            t0 = time()
            qtilphi = self.logPooling(phi) #vector with probability of each phi[i] belonging to qtilphi
            print "==> Done Running the Log Pooling (took %s seconds)\n"%(time()-t0)
            qtilphi = nan_to_num(qtilphi)
            #print 'max(qtilphi): ', max(qtilphi)
            if sum(qtilphi)==0:
                print 'Pooled prior on ouputs is null, please check your priors, and try again.'
                return 0
#      
        print "Calculating weights"
        po = Pool()
        jobs = [po.apply_async(fitfun, (phi[i],data)) for i in xrange(phi.shape[0])]
        w = [j.get() for j in jobs]
        po.close();po.join()
        w /=sum(w)
        w = 1-w
        
        w = nan_to_num(w)
        w = array(w)*qtilphi
        w /=sum(w)
        w = nan_to_num(w)
        print 'max(w): %s\nmin:%s\nmean(w): %s\nvar(w): %s'%(max(w), min(w), mean(w), var(w))
        if sum(w) == 0.0:
            print 'Resampling weights are all zero, please check your model or data.'
            return 0
        print "Resampling Thetas"
        t0 = time()
        j = 0
        while j < self.L: # Extract L samples from q1theta
            i=randint(0,w.size)# Random position of w and q1theta
            if random()<= w[i]:
                self.post_theta[j] = self.q1theta[i]# retain the sample according with resampling prob.
                j+=1
        print "==> Done Resampling (L=%s) priors (took %s seconds)"%(self.L,(time()-t0))

        self.done_running = True
        return 1

    def imp_sample(self,n,data, w):
        """
        Importance sampling
        
        :Returns:
            returns a sample of size n
        """
        #sanitizing weights
        print "Starting importance Sampling"
        w /=sum(w)
        w = nan_to_num(w)
        j=0
        k=0
        smp = copy.deepcopy(data[:n])
        while j < n:
            i = randint(0,w.size)# Random position of w
            if random() <= w[i]:
                smp[j] = data[j]
                j += 1
                
            k+=1
        print "Done imp samp."
        return smp
    
    def mcmc_run(self, data, t=1, likvariance=10,burnin=1000, nopool=False, method="MH" , constraints=[]):
        """
        MCMC based fitting

        :Parameters:
            - `data`: observed time series on the model's output
            - `t`: length of the observed time series
            - `likvariance`: variance of the Normal likelihood function
            - `nopool`: True if no priors on the outputs are available. Leads to faster calculations
            - `method`: Step method. defaults to Metropolis hastings
        """
        #self.phi = recarray((self.K,t),formats=['f8']*self.nphi, names = self.phi.dtype.names)
        ta = True if self.verbose==1 else False
        tc = True if self.verbose==1 else False
        if method == "MH":
            sampler = MCMC.Metropolis(self, self.K,self.K*10, data, t, self.theta_dists, self.q1theta.dtype.names, self.tlimits, like.Normal, likvariance, burnin, trace_acceptance=ta,  trace_convergence=tc, nchains=self.ntheta, constraints=[])
            sampler.step()
            self.phi = sampler.phi
            #self.mh(self.K,t,data,like.Normal,likvariance,burnin)
        elif method == 'dream':
            sampler = MCMC.Dream(self, self.K,self.K*10, data, t, self.theta_dists, self.q1theta.dtype.names, self.tlimits, like.Normal, likvariance, burnin, trace_acceptance=ta,  trace_convergence=tc, nchains=self.ntheta, constraints=[])
            sampler.step()
            self.phi = sampler.phi
        else:
            sys.exit("Invalid MCMC method!\nTry 'MH'.")
        self.done_running = 1
        return 1
        
    def _output_loglike(self, prop, data, likfun=like.Normal,likvar=1e-1, po=None):
        """
        Returns the log-likelihood of a simulated series
        
        :Parameters:
            - `prop`: Proposed output
            - `data`: Data against which proposal will be measured
            - `likfun`: Likelihood function
            - `likvar`: Variance of the likelihood function
            - `po`: Pool of processes for parallel execution
        
        :Types:
            - `prop`: array of shape (t,nphi) with series as columns.
            - `data`: Dictionary with keys being the names (as in phinames) of observed variables
            - `likfun`: Log likelihood function object
        """
        if isinstance(prop, numpy.recarray):
            prop= numpy.array(prop.tolist())
        t = prop.shape[0] #1 if model's output is a scalar, larger if it is a time series (or a set of them)
        lik=0
        for k in xrange(self.nphi): #loop on series
            if self.q2phi.dtype.names[k] not in data:
                continue#Only calculate liks of series for which we have data
            obs = array(data[self.q2phi.dtype.names[k]])
            if len(obs.shape)>1:#in case of more than one dataset
                obs = clearNaN(obs).mean(axis=1)
            if po != None:# Parallel version
                liks = [po.apply_async(likfun,(obs[p],prop[p][k],1./likvar)) for p in xrange(t) if not isnan(obs[p])]
                lik = nansum([l.get() for l in liks])
            else:
                liks = [likfun(obs[p],prop[p][k],1./likvar) for p in xrange(t)if not isnan(obs[p]) ]
                lik = nansum(liks)
#                if isnan(lik):
#                    pdb.set_trace()
        return lik
    
    def sir(self, data={}, t=1,variance=0.1, pool=False,savetemp=False):
        """
        Run the model output through the Sampling-Importance-Resampling algorithm.
        Returns 1 if successful or 0 if not.

        :Parameters:
            - `data`: observed time series on the model's output
            - `t`: length of the observed time series
            - `variance`: variance of the Normal likelihood function
            - `pool`: False if no priors on the outputs are available. Leads to faster calculations
            - `savetemp`: Boolean. create a temp file?
        """
        seed()
        phi = self.runModel(savetemp,t)
        # Do Log Pooling
        if not pool:
            qtilphi = ones(self.K)
        else:
            t0 = time()
            qtilphi = self.logPooling(phi) #vector with probability of each phi[i] belonging to qtilphi
            print "==> Done Running the Log Pooling (took %s seconds)\n"%(time()-t0)
            qtilphi = nan_to_num(qtilphi)
            print 'max(qtilphi): ', max(qtilphi)
            if sum(qtilphi)==0:
                print 'Pooled prior on ouputs is null, please check your priors, and try again.'
                return 0

#       Calculating the likelihood of each phi[i] considering the observed data
#        lik =self._output_loglike()
        lik = zeros(self.K)
        t0=time()
        po = Pool()
        for i in xrange(self.K):
            p = phi[i]
            l =self._output_loglike(p,data,  likvar=variance)
#            for n in data.keys():
#                if isinstance(data[n],list) and data[n] == []: 
#                    continue #no observations for this variable
#                elif isinstance(data[n],numpy.ndarray) and (not data[n].any()):
#                    continue #no observations for this variable
#                p = phi[n]
#                lik =self._output_loglike()
#                liklist=[po.apply_async(like.Normal,(data[n][m], j, 1./variance)) for m,j in enumerate(p[i])]
#                l=sum([p.get() for p in liklist])
            if i%self.K/10. ==0:
                print "Likelihood calculation progress: %s of %s done."%(i, self.K)
            lik[i]=l
        po.close()
        po.join()
        if self.viz:
            dtplot.clearFig();phiplot.clearFig();thplot.clearFig()
            dtplot.gp.xlabel('observed')
            dtplot.gp.ylabel('simulated')
            obs = [];sim =[]
            for n in data.keys():
                obs.append(data[n])
                sim.append(phi[n].mean(axis=0).tolist())
            dtplot.scatter(array(obs),array(sim),names=data.keys(),title='fit')
            phiplot.plotlines(array(sim),names=data.keys(),title='Model Output')
            thplot.plothist(self.q1theta, title='Input parameters',names=self.q1theta.dtype.names)
        print "==> Done Calculating Likelihoods (took %s seconds)"%(time()-t0)
        lr = nan_to_num(max(lik)/min(lik))
        print '==> Likelihood (min,mean,max,sum): ',min(lik),mean(lik),max(lik), sum(lik)
        print "==> Likelihood ratio of best run/worst run: %s"%(lr,)
#        Calculating the weights
        w = nan_to_num(qtilphi*lik)
        w = nan_to_num(w/sum(w))

        if not sum(w) == 0.0:
            j = 0
            t0 = time()
            maxw = 0;minw = max(w) #keep track of goodness of fit of phi
            while j < self.L: # Extract L samples from q1theta
                i=randint(0,w.size)# Random position of w and q1theta
                if random()*max(w)<= w[i]:
                    self.post_theta[j] = self.q1theta[i]# retain the sample according with resampling prob.
                    maxw = max(maxw,w[i])
                    minw = min(minw,w[i])
                    j+=1
                    if not j%100 and self.verbose:
                        print j, "of %s"%self.L
            self.done_running = True
            print "==> Done Resampling (L=%s) priors (took %s seconds)"%(self.L,(time()-t0))
            wr = maxw/minw
            print "==> Likelihood ratio of best/worst retained runs: %s"%(wr,)
            if wr == 1:
                print "==> Flat likelihood, trying again..."
                return 0
            print "==> Improvement: %s percent"%(100-100*wr/lr,)
        else:
            print 'Resampling weights are all zero, please check your model or data, and try again.\n'
            print '==> Likelihood (min,mean,max): ',min(lik),mean(lik),max(lik)
            print '==> RMS deviation of outputs: %s'%(basicfit(phi, data),)
            return 0
        return 1

        
    def runModel(self,savetemp,t=1, k=None):
        '''
        Handles running the model k times keeping a temporary savefile for
        resuming calculation in case of interruption.

        :Parameters:
            - `savetemp`: Boolean. create a temp file?
            - `t`: number of time steps
        
        :Returns:
            - self.phi: a record array of shape (k,t) with the results.
        '''
        if savetemp:
            CP.dump(self.q1theta,open('q1theta','w'))
        if not k:
            k = self.K
#       Running the model ==========================
        
                
        if os.path.exists('phi.temp'):
            self.phi,j = CP.load(open('phi.temp','r'))
        else:
            j=0
            self.phi = recarray((k,t),formats=['f8']*self.nphi, names = self.phi.dtype.names)
        def cb(r):
            '''
            callback function for the asynchronous model runs
            '''
            if t == 1:
                self.phi[r[1]] = (r[0][-1],)
            else:
                self.phi[r[1]] = [tuple(l) for l in r[0][-t:]]# #phi is the last t points in the simulation

        po = Pool()
        t0=time()
        for i in xrange(j,k):
            theta = [self.q1theta[n][i] for n in self.q1theta.dtype.names]
            r = po.apply_async(enumRun,(self.model,theta,i),callback=cb)
#            r = po.apply_async(self.model,theta)
#            if t == 1:
#                phi[i] = (r.get()[-1],)
#            else:
#                phi[i] = [tuple(l) for l in r.get()[-t:]]# #phi is the last t points in the simulation
            if i%100 == 0 and self.verbose:
                print "==> K = %s"%i
                if savetemp:
                    CP.dump((self.phi,i),open('phi.temp','w'))
        if savetemp: #If all replicates are done, clear temporary save files.
            os.unlink('phi.temp')
            os.unlink('q1theta')
        po.close()
        po.join()
        print "==> Done Running the K (%s) replicates (took %s seconds)\n"%(k,(time()-t0))
        
        return self.phi

def basicfit(s1,s2):
    '''
    Calculates a basic fitness calculation between a model-
    generated time series and a observed time series.
    it uses a Mean square error.

    :Parameters:
        - `s1`: model-generated time series. record array.
        - `s2`: observed time series. dictionary with keys matching names of s1

    :Return:
        Root mean square deviation between ´s1´ and ´s2´.
    '''
    if isinstance(s1, recarray):
#        print "==> is recarray!"
        assert isinstance(s2, dict)
        err = []
        for k in s2.keys():
            if k not in s1.dtype.names:
                continue
            ls1 = len(s1[k]) #handles the cases where data is slightly longer that simulated series.
#                print s2[k]
#                try:
#                pdb.set_trace()

            if len(s2[k].shape) >1:
                s2[k] = clearNaN(s2[k]).mean(axis=1)#nanmean(s2[k], axis=1)
            dif = s1[k]-s2[k][:ls1].astype(float)
            
            dif[isnan(dif)] = 0
            e = nanmean(dif**2)
#                except TypeError:
#                    print s1[k], s2[k]

            err.append(e) 
    elif isinstance(s1, list):
#        print "==> is List!"
        assert isinstance(s2, list) and len(s1) ==len(s2)
        s1 = array(s1).astype(float)
        s2 = array(s2).astype(float)
        if len(s2.shape) >1:
            s2 = clearNaN(s2).mean(axis=1)
        
        err = [nanmean((s-t)**2) for s, t in zip(s1, s2) ]
        #err = [sum((s-t)**2./t**2) for s, t in zip(s1, s2)]
    mse = nan_to_num(mean(err))
    return mse

def clearNaN(obs):
    """
    Loops through an array with data series as columns, and 
    Replaces NaNs with the mean of the other series.
    
    :Parameters:
        - `obs`: 2-dimensional numpy array
        
    :Returns:
        array of the same shape as obs
    """
    rows = obs.T.tolist() #tranpose to facilitate processing
    res  = zeros(obs.T.shape)
    for i, r in enumerate(rows):
        a = array(r)
        if not isnan(a).any():
            continue
        else:
            for j, e in enumerate(r):
                if isnan(e):
                    c = obs.T[:, j].tolist()
                    c.pop(i)
                    r[j] = nanmean(c)
                    r = nan_to_num(r)
            res[i, :] = r
                
    return res.T

def enumRun(model,theta,k):
    """
    Returns model results plus run number.

    :Parameters:
        - `model`: model callable
        - `theta`: model input list
        - `k`: run number
        
    :Return:
        - res: result list
        - `k`: run number
    """
    res =model(theta)
    return (res,k)

def model_as_ra(theta, model, phinames):
    """
    Does a single run of self.model and returns the results as a record array
    """
    theta = list(theta)
    nphi = len(phinames)
    r = model(theta)
    res = recarray(r.shape[0],formats=['f8']*nphi, names = phinames)
    for i, n in enumerate(res.dtype.names):
        res[n] = r[:, i]
    return res

def model(theta, n=1):
    """
    Model (r,p0, n=1)
    Simulates the Population dynamic Model (PDM) Pt = rP0
    for n time steps.
    P0 is the initial population size. 
    Example model for testing purposes.
    """
    r, p0 = theta
#    print "oi"
    Pt = zeros((n, 1), float) # initialize the output vector
    P = p0
    for i in xrange(n):
        Pt[i] = r*P
        P = Pt[i]
    
    return Pt


def plotRaHist(arr):
    '''
    Plots a record array
    as a panel of histograms
    '''
    nv = len(arr.dtype.names)
    fs = (numpy.ceil(numpy.sqrt(nv)),numpy.floor(numpy.sqrt(nv))+1) #figure size
    P.figure()
    for i,n in enumerate(arr.dtype.names):
        P.subplot(floor(sqrt(nv)),floor(sqrt(nv))+1,i+1)
        P.hist(arr[n],bins=50, normed=1, label=n)
        P.legend()


def main2():
    start = time()
    #Me = Meld(K=5000,L=1000,model=model, ntheta=2,nphi=1,verbose=False,viz=False)
    Me.setTheta(['r','p0'],[stats.uniform,stats.uniform],[(2,4),(0,5)], [(0, 10), (0, 10)])
    Me.setPhi(['p'],[stats.uniform],[(6,9)],[(0,10)])
    #Me.add_data(normal(7.5,1,400),'normal',(6,9))
    #Me.run()
    Me.sir(data ={'p':[7.5]} )
    pt,pp = Me.getPosteriors(1)
    end = time()
    print end-start, ' seconds'
    plotRaHist(pt)
    plotRaHist(pp)
    P.show()
    
    
def mh_test():
    start = time()
    #Me = Meld(K=5000,L=1000,model=model, ntheta=2,nphi=1,verbose=False,viz=False)
    Me.setTheta(['r','p0'],[stats.uniform,stats.uniform],[(2,4),(0,5)], [(0, 10), (0, 10)])
    Me.setPhi(['p'],[stats.uniform],[(6,9)],[(0,10)])
    Me.mcmc_run(data ={'p':[7.5]},burnin=1000)
    pt,pp = Me.getPosteriors(1)
    end = time()
    print end-start, ' seconds'
    plotRaHist(pt)
    plotRaHist(pp)
    P.show()
    

if __name__ == '__main__':
    Me = Meld(K=5000,L=1000,model=model, ntheta=2,nphi=1,verbose=True,viz=False)
    #mh_test()
    main2()
     

