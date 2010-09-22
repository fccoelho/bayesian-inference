# To change this template, choose Tools | Templates
# and open the template in the editor.
"""
Module implementing MCMC samplers 

    - Metropolis: Adaptive Metropolis Hastings sampler
    - Dream: DiffeRential Evolution Adaptive Markov chain sampler
"""
import sys
import time
import pdb
import cython
import xmlrpclib
import logging
from multiprocessing import Pool
from random import sample

import numpy as np
from liveplots.xmlrpcserver import rpc_plot
from numpy import array, mean,isnan,  nan_to_num, var, sqrt, inf, exp, greater, less, identity, ones, zeros, floor, log, recarray, nan
from numpy.random import random,  multivariate_normal,  multinomial,  rand
from scipy.stats import cov,  uniform, norm, scoreatpercentile


__author__="fccoelho"
__date__ ="$09/12/2009 10:44:11$"
__docformat__ = "restructuredtext en"

logger = logging.getLogger('BIP.MCMC')

def timeit(method):
    """
    Decorator to time methods
    """
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print '%r  %2.2f sec' % \
              (method.__name__ , te-ts)
        return result

    return timed



class _Sampler(object): 
    '''
    Base classe for all samplers
    Holds common logic and 
    '''
    _po = None
    _dimensions = None #cache for dimensions
    trace_acceptance = False
    trace_convergence = False
    seqhist = None
    liklist = []
    e = 1e-20 #very small number used in the proposal covariance function calculation
    _j=-1
    _R = np.inf #Gelman Rubin Convergence
    def __init__(self,  parpriors=[],  parnames = []):
        self.parpriors = parpriors
        self.parnames = parnames

    @property
    def best_prop_index(self):
        '''
        Returns the index of the best fitting proposal, i.e.,
        the one which with max Likelihood
        '''
        if not self.liklist:
            return 0
        return self.liklist.index(max(self.liklist))

    @property
    def DIC(self):
        """
        Calculates  the deviance information criterion
        """
        D = -2*array(self.liklist)
        Dbar = nan_to_num(D).mean()
        meanprop = array([self.meld.post_phi[i].mean(axis=0) for i in self.meld.post_phi.dtype.names])
        pd = Dbar+2*self.meld._output_loglike(meanprop.T, self.data, self.likfun, self.likvariance)
        DIC = pd +Dbar
        return DIC

    @property
    def dimensions(self):
        if not self._dimensions:
            self._dimensions = len(self.parpriors)
        return self._dimensions

    @property
    def po(self):
        '''
        Pool of processes for parallel execution of tasks
        Remember to call self.term_pool() when done.
        '''
        if self._po == None:
            self._po = Pool()
        else:
            if self._po._state:
                self._po = Pool() #Pool has been terminated
        return self._po

    def shut_down(self,reason=''):
        '''
        Finalizes the sampler, nicely closing the resources allocated

        :Parameters:
            - `reason`: comment stating why the sampling is being shutdown.
        '''
        try:
            self.term_pool()
        except OSError:
            pass
        self.pserver.close_plot()
        self.pserver2.close_plot()
        if reason:
            logger.info(reason)


    def term_pool(self):
        if self._po == None:
            return
        if not self._po._state: #Pool needs terminating
            self._po.close()
            self._po.join()
            self._po = None

    def gr_R(self, end, start):
        if self._j == end:
            return self._R
        else:
            self.gr_convergence(end, start)
            self._j = end
            return self._R

    def gr_convergence(self, relevantHistoryEnd, relevantHistoryStart):
        """
        Gelman-Rubin Convergence
        """
        start = relevantHistoryStart
        end = relevantHistoryEnd
        N = end - start
        if N==0:
            self._R = np.inf*np.ones(self.nchains)
            return
        N = min(min([len(self.seqhist[c]) for c in range(self.nchains)]), N)
        seq = [self.seqhist[c][-N:] for c in range(self.nchains)]
        sequences = array(seq) #this becomes an array (nchains,samples,dimensions)
        variances  = var(sequences,axis = 1)#array(nchains,dim)
        means = mean(sequences, axis = 1)#array(nchains,dim)
        withinChainVariances = mean(variances, axis = 0)
        betweenChainVariances = var(means, axis = 0) * N
        varEstimate = (1 - 1.0/N) * withinChainVariances + (1.0/N) * betweenChainVariances
        self._R = sqrt(varEstimate/ withinChainVariances)

    @np.vectorize
    def _accept(self, last_lik,  lik):
        """
        Decides whether to accept a proposal
        """
        if last_lik == None: last_lik = -inf
        # liks are logliks
        if lik == -inf:#0:
            return 0
        if last_lik >-inf:#0:
            alpha = min( exp(lik-last_lik), 1)
            #alpha = min(lik-last_lik, 1)
        elif last_lik == -inf:#0:
            alpha = 1
        else:
            return 0
            raise ValueError("Negative likelihood!?!")
#        print "last_lik, lik, alpha: ",  last_lik, lik, alpha
        if random() < alpha:
            return 1
        else:
            return 0

    def setup_xmlrpc_plotserver(self):
        """
        Sets up the server for real-time chain watch
        """
        p=0;p2=0
        while p==0 or p2 == 0:
            p = rpc_plot()
            p2 = rpc_plot(hold=1)
        self.pserver = xmlrpclib.ServerProxy('http://localhost:%s'%p)
        self.pserver2 = xmlrpclib.ServerProxy('http://localhost:%s'%p2)
        
    def shutdown_xmlrpc_plotserver(self):
        self.pserver.flush_queue()
        self.pserver.shutdown()
        self.pserver2.flush_queue()
        self.pserver2.shutdown()
        
    def _every_plot(self):
        """
        plotting function for generating a plot at every step
        """
        pass

    def _watch_chain(self, j):
        if j<100:
            return
        self.gr_convergence(j, j-100)
        print "Gelman-Rubin's R: ", self._R
        self.pserver.clearFig()
        thin = j//500 if j//500 !=0 else 1 #history is thinned to show at most 500 points, equally spaced over the entire range
        chaindata = self.history[:j:thin].T.tolist()
        obs = [];lbs = []
        for k, d in self.data.items():
            if len(d.shape)>1:
                obs += [nan_to_num(i).tolist() for i in d.T]
                lbs += [k+str(i) for i in range(d.shape[1])]
            else:
                obs += nan_to_num(d).tolist()
                lbs += [k]
        self.pserver.lines(chaindata,range(j-(len(chaindata[0])), j), self.parnames, "Chain Progress.",'points' , 1)
        self.pserver2.lines(obs,[],lbs, "Fit", 'points' )
        s = j-50 if 3*j//4<50 else 3*j//4
        #series = [self.phi[k][s:j].mean(axis=0).tolist() for k in self.data.keys()]
        series = [mean(self.phi[k][s:j], axis=0).tolist() for k in self.data.keys()]
        bi = self.liklist.index(max(self.liklist[s:j])) #index of the best fit
        series = [mean(self.phi[k][s:j], axis=0).tolist() for k in self.data.keys()]
        best = [self.phi[k][bi].tolist() for k in self.data.keys()]
        mlabels = ['Mean '+l for l in self.data.keys()]
        blabels = ['Best '+l for l in self.data.keys()]
        self.pserver2.lines(series,[],mlabels, "Mean fit of last %s samples"%(j-s), 'lines' )
        self.pserver2.lines(best,[],blabels, "mean and best fit of last %s samples"%(j-s), 'lines' )
        self.pserver2.clearFig()
        #TODO: Implement plot of best fit simulation against data

    def _tune_likvar(self, ar):
        try:
            self.arhist.append(ar)
        except AttributeError:
            self.tsig = 1
            self.tstep = .05
            self.arhist = [ar]
        dev = (0.35-ar)**2
        if dev > 0.02:
            self.likvariance *= 1+self.tsig *(.5*(np.tanh(8*dev-3)+1))
        else: return #ar at target, don't change anything
        improv = (0.35-mean(self.arhist[-5:-1]))**2 - (0.35-ar)**2
        if improv < 0:
            self.tsig *= -1 #change signal if AR is not improving
            self.tstep = .05 #reset to small steps if changing direction
        elif improv  > 0 and improv <.01:
            if random() <.05: #1 in 20 chance to change direction if no improvements
                self.tsig *= -1 #change signal if AR is not improving
        elif improv > 0.01:
            self.tstep *= 0.97 #reduce step if approacching sweet spot

#    @np.vectorize
    def check_constraints(self, theta):
        """
        Check if given theta vector complies with all constraints
        
        :Parameters:
            - `theta`: parameter vector
            
        :Returns:
            True if theta passes all constraints, False otherwise
        """
        if not self.constraints:
            return True
        r = array([c(theta) for c in self.constraints])
        return r.all()
    
    def _propose(self, step, po=None):
        """
        Generates proposals.
        returns two lists

        :Parameters:
            - `step`: Position in the markov chain history.
            - `po`: Process pool for parallel proposal generation

        :Returns:
            - `theta`: List of proposed self.dimensional points in parameter space
            - `prop`: List of self.nchains proposed phis.
        """
        thetalist = []
        proplist = []
        initcov = np.identity(self.dimensions)
        for c in range(self.nchains):
            if step <= 1 or self.seqhist[c] ==[]: 
                #sample from the priors
                while 1:
                    theta = [self.parpriors[dist]() for dist in self.parnames]
                    if not self.check_constraints(theta):
                        continue
                    if sum ([int(greater(t, self.parlimits[i][0]) and less(t, self.parlimits[i][1])) for i, t in enumerate(theta)]) == self.dimensions:
                        break
                self.lastcv = initcov #assume no covariance at the beginning
            else:
                #use gaussian proposal
                if step%10==0 and len(self.seqhist[c]) >=10: #recalculate covariance matrix only every ten steps
                    cv = self.scaling_factor*cov(array(self.seqhist[c][-10:]))+self.scaling_factor*self.e*identity(self.dimensions)
                    self.lastcv = cv
                else:
                    cv = self.lastcv
                while 1:
                    theta = multivariate_normal(self.seqhist[c][-1],cv, size=1).tolist()[0]
                    if sum ([int(greater(t, self.parlimits[i][0]) and less(t, self.parlimits[i][1])) for i, t in enumerate(theta)]) == self.dimensions:
                        break
            thetalist.append(theta)
        if po:
            proplis = [po.apply_async(model_as_ra, (t, self.meld.model, self.meld.phi.dtype.names)) for t in thetalist]
            proplist = [job.get() for job in proplis]
        else:
            proplist = [model_as_ra(t, self.meld.model, self.meld.phi.dtype.names) for t in thetalist]
        propl = [p[:self.t] for p in proplist]
        return thetalist,propl

#TODO: remove dependency on the meld object
class Metropolis(_Sampler):
    """
    Standard random-walk Metropolis Hastings sampler class
    """
    def __init__(self, meldobj, samples, sampmax,  data, t, parpriors, parnames, parlimits, likfun, likvariance, burnin, **kwargs):
        """
        MCMC based fitting

        :Parameters:
            - `samples`: Number of samples to obtain
            - `sampmax`: Maximum number of samples drawn.
            - `data`: observed time series on the model's output
            - `t`: length of the observed time series
            - `parpriors`: Dictionary with frozen distributions objects as values and parnames as keys
            - `parnames`: List of parameter names
            - `parlimits`: list of tuples with (min,max) for every parameter.
            - `likfun`: Likelihood function
            - `likvariance`: variance of the Normal likelihood function
            - `burnin`: Number of burnin samples
        """
        self.salt_band = 0.05
        self.samples = samples
        self.sampmax = sampmax
        self.parpriors = parpriors
        self.parnames = parnames
        self.parlimits = parlimits
        self.likfun = likfun
        self.likvariance = likvariance
        self.data = data
        self.meld = meldobj
        self.t = t
        self.burnin = burnin
        self.nchains = 1 
        self.phi = np.recarray((self.samples+self.burnin,t),formats=['f8']*self.meld.nphi, names = self.meld.phi.dtype.names)
        self.scaling_factor = (2.38**2)/self.dimensions
        self.e = 1e-20
        if kwargs:
            for k, v in kwargs.iteritems():
                exec('self.%s = %s'%(k, v))
        self.nchains = 1 
        # Combined history of accepted samples
        self.history = np.zeros((self.nchains*(samples+self.burnin), self.dimensions)) 
        #complete history of all chains as a dictionary with keys as integer ids of the chains
        self.seqhist = dict([(i, [])for i in range(self.nchains)])
        #self.seqhist = np.zeros((self.nchains, self.dimensions, samples+self.burnin))
        self.setup_xmlrpc_plotserver()

    def _propose(self, step, po=None):
        """
        Generates proposals.
        returns two lists

        :Parameters:
            - `step`: Position in the markov chain history.
            - `po`: Process pool for parallel proposal generation

        :Returns:
            - `theta`: List of proposed self.dimensional points in parameter space
            - `prop`: List of self.nchains proposed phis.
        """
        po=None
        thetalist = []
        proplist = []
        initcov = identity(self.dimensions)
        if self.meld.initheta and step <= 1:
            #start from user-defined point in parameter space.
            for i in range(self.nchains):
                thetalist.append(self.meld.initheta)
            self.lastcv = initcov #assume no covariance at the beginning
        else:
            for c in range(self.nchains):
                off = 0
                if step <= 1 or self.seqhist[c] ==[]: 
                    #sample from the priors
                    while off<50:
                        theta = [self.parpriors[par].rvs() for par in self.parnames]
                        if not self.check_constraints(theta):
                            continue
                        if sum ([int(t>= self.parlimits[i][0] and t<= self.parlimits[i][1]) for i, t in enumerate(theta)]) == self.dimensions:
                            break
                        off+=1
                    if off ==50:#try a compromising proposal
                        theta = self.seqhist[c][-1] #last accepted proposal for this chain
    #                print "off:" , off
                    self.lastcv = initcov #assume no covariance at the beginning
                else:
                    #use gaussian proposal
                    if step%10==0 and len(self.seqhist[c]) >=10: #recalculate covariance matrix only every ten steps
                        cv = self.scaling_factor*cov(array(self.seqhist[c][-10:]))+self.scaling_factor*self.e*identity(self.dimensions)
                        self.lastcv = cv
                    else:
                        cv = self.lastcv
                    #print self.parlimits
                    while off<50:
                        theta = multivariate_normal(self.seqhist[c][-1],cv, size=1).tolist()[0]
                        if sum ([int(t>= self.parlimits[i][0] and t<= self.parlimits[i][1]) for i, t in enumerate(theta)]) == self.dimensions:
                            break
                        off+=1
                    if off ==50: #try a compromising proposal
                        theta = self.seqhist[c][-1] #last accepted proposal for this chain
                    #print "off:" , off
                thetalist.append(theta)
        if po:
            proplis = [po.apply_async(model_as_ra, (t, self.meld.model, self.meld.phi.dtype.names)) for t in thetalist]
            proplist = [job.get() for job in proplis]
        else:
            proplist = [model_as_ra(t, self.meld.model, self.meld.phi.dtype.names) for t in thetalist]
        propl = [p[:self.t] for p in proplist]
        return thetalist,propl

    def step(self,  nchains=1):
        """
        Does the actual sampling loop.
        """
        ptheta = recarray(self.samples+self.burnin,formats=['f8']*self.dimensions, names = self.parnames)
        i=0;j=0;rej=0;ar=0 #total samples,accepted samples, rejected proposals, acceptance rate
        last_lik = None
        while j < self.samples+self.burnin:
            print j
            self.meld.current_step = j
            if self.meld.stop_now:
                return self.shut_down('user interrupted')
            #generate proposals
            theta,prop = self._propose(j, self.po)
            #calculate likelihoods
            lik = [self.meld._output_loglike(p, self.data, self.likfun, self.likvariance, self.po) for p in prop]

#            print "lik:" , lik,  last_lik,  j
            accepted = self._accept(self, last_lik, lik)# have to include self in the call because method is vectorized.
#            print "acc:", accepted,  theta
            #Decide whether to accept proposal
            if last_lik == None: #on first sample
                last_lik = lik
                continue
            i +=self.nchains
            if sum(accepted) < self.nchains:
                rej += self.nchains-sum(accepted) #adjust rejection counter
                if i%100 == 0: 
                    ar = (i-rej)/float(i)
                    self._tune_likvar(ar)
                    if self.trace_acceptance:
                        print "--> %s: Acc. ratio: %s"%(rej, ar)
            # Store accepted values
#            print "nchains:", self.nchains
            for c, t, pr,  a in zip(range(self.nchains), theta, prop, accepted): #Iterates over the results of each chain
                #if not accepted repeat last value
                if not a:
                    continue
                self.history[j, :] = t 
                self.seqhist[c].append(t)
                #self.seqhist[c, :, j] = t 
                self.phi[j] = pr[0] if self.t==1 else [tuple(point) for point in pr]
                ptheta[j] = tuple(t)
                self.liklist.append(lik[c])
                if j == self.samples+self.burnin:break
                j += 1 #update accepted sample counter 
            #print j,  len(self.seqhist[0])
            if j%100==0 and j>0:
                if self.trace_acceptance:
                    print "++>%s,%s: Acc. ratio: %s"%(j,i, ar)
                    self._watch_chain(j)
                if self.trace_convergence: print "++> %s: Likvar: %s\nML:%s"%(j, self.likvariance, np.max(self.liklist) )
#            print "%s\r"%j
            last_lik = lik
            last_prop = prop
            last_theta = theta
            ar = (i-rej)/float(i)
            if self.meld.verbose ==2 and j>10:
                self.meld.current_plot(self.phi, self.data, self.best_prop_index, step=j)
        self.term_pool()
        self.meld.post_theta = ptheta[self.burnin:]
        self.meld.post_phi = self.phi[self.burnin:]
        self.meld.post_theta = ptheta#self._imp_sample(self.meld.L,ptheta,liklist)
        self.meld.likmax = max(self.liklist)
        self.meld.DIC = self.DIC
        print "Total steps(i): ",i,"rej:",rej, "j:",j
        print ">>> Acceptance rate: %s"%ar
        self.shut_down('Finished normally')
        return 1


    def _rms_fit(self, s1, s2):
        '''
        Calculates a basic fitness calculation between a model-
        generated time series and a observed time series.
        It uses a normalized RMS variation.

        :Parameters:
            - `s1`: model-generated time series. 
            - `s2`: observed time series. dictionary with keys matching names of s1
        :Types:
            - `s1`: Record array or list.
            - `s2`: Dictionary or list

        s1 and s2 can also be both lists of lists or lists of arrays of the same length.

        :Return:
            Inverse of the Root mean square deviation between `s1` and `s2`.
        '''
        if isinstance(s1, np.recarray):
            assert isinstance(s2, dict)
            err = []
            for k in s2.keys():
                e = np.sqrt(np.mean((s1[k]-s2[k])**2.))
                err.append(e) 
        if isinstance(s1, list):
            assert isinstance(s2, list) and len(s1) ==len(s2)
            err = [np.sqrt(np.mean((s-t)**2.)) for s, t in zip(s1, s2)]
        rmsd = np.mean(err)
        fit = 1./rmsd #fitness measure
#        print "rmsd, fit, err: ", rmsd,fit, err
        if fit ==np.inf:
            sys.exit()
        return fit #mean r-squared

    def _imp_sample(self,n,data, w):
        """
        Importance sampling

        :Parameters:
            - `n`: Number of samples to return
            - `data`: record array (containing on or more vectors of data) to be resampled
            - `w`: Weight vector
        :Returns:
            returns a sample of size n
        """
        #sanitizing weights
        print "Starting importance Sampling"
        w /= sum(w)
        w = np.nan_to_num(w)
        j=0
        k=0
        nvar = len(data.dtype.names)
        smp = np.recarray(n,formats = [data.dtype.descr[0][1]]*nvar,names = data.dtype.names)
        #smp = copy.deepcopy(data[:n])
        while j < n:
            i = np.random.randint(0,w.size)# Random position of w
            if np.random.random() <= w[i]:
                smp[j] = data[j]
                j += 1
            k+=1
        print "Done importance sampling."
        return smp

#    def _watch_chain(self, j):
#        if j<100:
#            return
#        self.gr_convergence(j, j-100)
#        self.pserver.clearFig()
#        thin = j//500 if j//500 !=0 else 1 #history is thinned to show at most 500 points, equally spaced over the entire range
#        data = self.history[:j:thin].T.tolist()
#        self.pserver.plotlines(data,range(j-(len(data[0])), j), self.parnames, "Chain Progress. GR Convergence: %s"%self._R,'points' , 1) 

    def _add_salt(self,dataset,band):
        """
        Adds a few extra uniformly distributed data 
        points beyond the dataset range.
        This is done by adding from a uniform dist.

        :Parameters:
            - `dataset`: vector of data
            - `band`: Fraction of range to extend: [0,1[
        :Returns:
            Salted dataset.
        """
        dmax = max(dataset)
        dmin = min(dataset)
        drange = dmax-dmin
        hb = drange*band/2.
        d = numpy.concatenate((dataset,stats.uniform(dmin-hb,dmax-dmin+hb).rvs(self.K*.05)))
        return d


def model_as_ra(theta, model, phinames):
    """
    Does a single run of self.model and returns the results as a record array
    """
    theta = list(theta)
    nphi = len(phinames)
    r = model(theta)
    res = np.recarray(r.shape[0],formats=['f8']*nphi, names = phinames)
    for i, n in enumerate(res.dtype.names):
        res[n] = r[:, i]
    return res

class Dream(_Sampler):
    '''
    DiffeRential Evolution Adaptive Markov chain sampler
    '''
    def __init__(self, meldobj, samples, sampmax, data, t , parpriors, parnames, parlimits,likfun, likvariance, burnin, thin = 5, convergenceCriteria = 1.1,  nCR = 3, DEpairs = 1, adaptationRate = .65, eps = 5e-6, mConvergence = False, mAccept = False, **kwargs):
        self.meld = meldobj
        self.samples = samples
        self.sampmax = sampmax
        self.data = data
        self.t = t
        self.parpriors = parpriors
        self.parnames = parnames
        self.parlimits = parlimits
        self.likfun = likfun
        self.likvariance = likvariance
        self.burnin = burnin
        self.nchains = len(parpriors) 
        self.phi = np.recarray((self.samples+self.burnin,t),formats=['f8']*self.meld.nphi, names = self.meld.phi.dtype.names)
        self.nCR = nCR
        self.DEpairs = DEpairs
        self.delayRej = 1
        if kwargs:
            for k, v in kwargs.iteritems():
                exec('self.%s = %s'%(k, v))
        self._R = array([2]*self.nchains) #initializing _R
        self.maxChainDraws = floor(samples/self.nchains)
        #initialize the history arrays
        # History of log posterior probs for all chains
        self.omega = zeros((self.samples+self.burnin, self.nchains))
        # Combined history of accepted samples
        self.history = zeros((self.nchains*(samples+self.burnin), self.dimensions))
        self.seqhist = dict([(i, [])for i in range(self.nchains)])
        #self.sequenceHistories = np.zeros((self.nchains, self.dimensions, self.maxChainDraws))
        # initialize the temporary storage vectors
        self.currentVectors = zeros((self.nchains, self.dimensions))
        self.currentLiks = ones(self.nchains)*-inf
        self.scaling_factor = 2.38/sqrt(2*DEpairs*self.dimensions)
        self.setup_xmlrpc_plotserver()

    def _det_outlier_chains(self, step):
        """
        Determine which chains are outliers
        """
        means = self.omega[step//2:step,:].mean(axis=0)
        q1 = scoreatpercentile(self.omega[step//2:step,:], 25)
        q3 = scoreatpercentile(self.omega[step//2:step,:], 75)
        iqr = q3-q1
        outl = means<q1-2*iqr
        return outl
#    @timeit
    def delayed_rejection(self, xi, zi, pxi, zprob):
        """
        Generates a second proposal based on rejected proposal xi
        """
        k=.3 #Deflation factor for the second proposal
        cv = self.scaling_factor*cov(xi)+self.scaling_factor*self.e*identity(self.dimensions)
        o=0
        while o<50:
            zdr = multivariate_normal(xi,k*cv,1).tolist()[0]
            if not self.check_constraints(zdr): continue
            if sum ([t>= self.parlimits[i][0] and t <= self.parlimits[i][1] for i, t in enumerate(zdr)]) == self.dimensions:
                break
            o+=1
        if not sum ([t>= self.parlimits[i][0] and t <= self.parlimits[i][1] for i, t in enumerate(zdr)]) == self.dimensions:
            return xi, 0, 0, 0, 0
        propphi_zdr = self._prop_phi([zdr])
#        print propphi_zdr, zdr
        zdrprob,  zdrlik = self._get_post_prob([zdr],propphi_zdr)
        alpha2 = min(zdrprob[0]*(1-self._alpha1(self,zdrprob[0],zprob))/pxi*(1-self._alpha1(self, pxi, zprob)), 1)
        acc = 0; lik = 0; pr = 0; prop = 0
        if random()< alpha2:
            xi = zdr
            acc = 1
            liks = zdrlik
            pr = zdrprob[0]
            prop = propphi_zdr
        return xi, acc, lik, pr, prop

    @np.vectorize
    def _alpha1(self, p1, p2):
        """
        Returns the Metropolis acceptance probability:
        alpha1(p1,p1) = min(1,p1/p2) if p2 >-np.inf else 1

        :Parameters:
            - `p1`: log probability
            - `p2`: log probability
        """
        if p2 == None: p2 = -inf
        # ps are log probabilities
        if p2 >-inf:#np.exp(p2)>0
            alpha = min( exp(p1-p2), 1)
        elif p2 == -inf:#np.exp(p2)==0
            alpha = 1
        else:
            print "proposal's logP: ", p2
            alpha = 0
        return alpha
        

    def update_CR_dist(self):
        t = 1
        Lm = 0 
        pm =1./self.nCR

        for i in range(self.nchains):
            m = multinomial(1, [pm]*self.nCR).nonzero()[0][0]+1
            CR = float(m)/self.nCR
            Lm +=1
            #TODO: finish implementing this

    def _prop_initial_theta(self, step):
        """
        Generate Theta proposals from priors
        """
        if self.meld.initheta:
            #start from user-defined point in parameter space.
            return [self.meld.initheta for i in range(self.nchains)]
            
        thetalist = []
        initcov = identity(self.dimensions)
        for c in range(self.nchains):
            #sample from the priors
#            while 1:
            theta = array([self.parpriors[par].stats(moments='m') for par in self.parnames])
#                if sum ([int(t>= self.parlimits[i][0] and t<= self.parlimits[i][1]) for i, t in enumerate(theta)]) == self.dimensions:
#                    break
            self.lastcv = initcov #assume no covariance at the beginning

            thetalist.append(theta.tolist())
        return thetalist 
        
#    @timeit
    def _prop_phi(self, thetalist, po=None):
        """
        Returns proposed Phi derived from theta
        """
        if po:
            propl = [po.apply_async(model_as_ra, (t, self.meld.model, self.meld.phi.dtype.names)) for t in thetalist]
            proplist = [job.get()[:self.t]  for job in propl]
        else:
            proplist = [model_as_ra(t, self.meld.model, self.meld.phi.dtype.names)[:self.t] for t in thetalist]
        return proplist
        
#    @timeit
    def _chain_evolution(self, proptheta,  propphi, pps, liks):
        """
        Chain evolution as describe in ter Braak's Dream algorithm.
        """
        CR = 1./self.nCR
        b = [(l[1]-l[0])/10. for l in self.parlimits]
        delta = (self.nchains-1)//2 if self.nchains >2 else 1
        gam = 2.38/sqrt(2*delta*self.dimensions)
        zis = []
        for c in xrange(self.nchains):
            o = 0
            while 1: #check constraints
                e = [uniform(-i, 2*i).rvs() for i in b]
                eps = [norm(0, i).rvs() for i in b]
                others = [x for i, x in enumerate(proptheta) if i !=c]
                dif = zeros(self.dimensions)
                for d in range(delta):
                    d1, d2 = sample(others, 2)
                    dif+=array(d1)-array(d2)
                zi = array(proptheta[c])+(ones(self.dimensions)+e)*gam*dif+eps
                #revert offlimits proposals
                for i in xrange(len(zi)):
                    if zi[i]<= self.parlimits[i][0] or zi[i]>= self.parlimits[i][1]:# or isnan(zi):
                        zi[i] = proptheta[c][i]
                #Cross over
                for i in xrange(len(zi)): 
                    zi[i] = proptheta[c][i] if rand() < 1-CR else zi[i]
                zis.append(zi)
                if self.check_constraints(zi): 
                    break
            
        #get the associated Phi's
        if isnan(zis).any():
            pdb.set_trace()
        propphi_z = self._prop_phi(zis, self.po)
        zprobs,  zliks = self._get_post_prob(zis, propphi_z)
        prop_evo = [0]*self.dimensions
        liks_evo = [0]*self.dimensions
        
        evolved = [0]*self.dimensions #evolved Theta
        prop_evo = [0]*self.nchains
        liks_evo = [0]*self.dimensions
        pps_evo = zeros(self.nchains) #posterior probabilities
        accepted = self._accept(self, pps, zprobs)#have to pass self because method is vectorized
        
        # Do Delayed rejection with the chains that got rejected
        # and store results.
        i = 0
        for z, x in zip(zis, proptheta):
            if accepted[i]:
                evolved[i] = z
                prop_evo[i] = propphi_z[i]
                pps_evo[i] = zprobs[i]
                liks_evo[i] = zliks[i]
                #self.liklist.append(zliks[i])
            else:
                th2,acc,lk,pr,prop = self.delayed_rejection(x,z,pps[i],zprobs[i])
                if acc:
                    accepted[i] = 1
                    #self.liklist.append(lk)
                evolved[i] = th2
                prop_evo[i] = prop if acc else propphi[i]
                liks_evo[i] = lk if acc else liks[i]
                try:
                    pps_evo[i] = pr if acc else pps[i]
                except TypeError: #when pps == None
                    pps_evo[i] = -inf
            i += 1
        return evolved, prop_evo, pps_evo, liks_evo, accepted
        
#    @timeit
    def _get_post_prob(self, theta, prop, po = None):
        '''
        Calculates the posterior probability for the proposal of each chain

        :Parameters:
            - `theta`: list of nchains thetas
            - `prop`: list of nchains phis
            - `po`: Pool of processes

        :Returns:
            - `posts`: list of log posterior probabilities of length self.nchains
            - `listoliks`: list of log-likelihoods of length self.nchains
        '''
        pri = 1
        pris = []
        for c in xrange(len(theta)):#iterate over chains
            for i in xrange(len(theta[c])): #iterate  over parameters
                try:
                    pri *= self.parpriors[self.parnames[i]].pdf(theta[c][i])
                except AttributeError: #in case distribution is discrete
                    pri *= self.parpriors[self.parnames[i]].pmf(theta[c][i])
            pris.append(pri)
        if po:
            listol = [po.apply_async(self.meld._output_loglike, (p, self.data, self.likfun, self.likvariance)) for p in prop]
            listoliks = [l.get() for l in listol]
            self.term_pool()
        else:
            listoliks = [self.meld._output_loglike(p, self.data, self.likfun, self.likvariance) for p in prop]
#        Multiply by prior values to obtain posterior probs
#        Actually sum the logs
        posts = (log(array(pris))+array(listoliks)).tolist()
        
        if isnan(posts).any():
            print "\nLikelihoods returned some NaNs. Dropping to debug mode:\n"
            pdb.set_trace()
        return posts, listoliks
    
    def step(self):
        """
        Does the actual sampling loop.
        """
        ptheta = recarray(self.samples+self.burnin,formats=['f8']*self.dimensions, names = self.parnames)
        i = 0;j=0;rej=0;ar=0 #total samples,accepted samples, rejected proposals, acceptance rate
        last_pps = None
        t0=time.time()
        while j < self.samples+self.burnin:
            self.meld.current_step = j
            if self.meld.stop_now:
                return self.shut_down('user interrupted')
            #generate proposals
            if j == 0:
                theta = self._prop_initial_theta(j)
                prop = self._prop_phi(theta, self.po)
                pps, liks = self._get_post_prob(theta, prop)
            else:
                theta = [self.seqhist[c][-1] for c in range(self.nchains)]
                prop = self._prop_phi(theta, self.po)
                #pps = last_pps
                #liks = last_liks
            # Evolve chains
#            while sum(self._R <=1.2)<self.nchains:
            theta, prop, pps, liks, accepted = self._chain_evolution(theta, prop, pps, liks)
            #storing log post probs
#            print self.omega.shape,  pps.shape
            self.omega[j, :] = pps
            #Compute GR R
            self.gr_R(j, -j//2)
#            if sum(self._R <=1.2)==self.nchains:
#                print "Converged on all dimensions"
#                print j, self._R
            #Update last_lik
            if last_pps == None: #on first sample
                last_pps = pps
                #last_liks = liks
                continue
            i +=self.nchains
            if sum(accepted) < self.nchains:
                ar = (i-rej)/float(i)
                rej += self.nchains-sum(accepted) #adjust rejection counter
#                print "==> Acc. ratio: %2.2f"%ar
                if i%100 == 0: 
                    self._tune_likvar(ar)
                    if self.trace_acceptance:
                        print "--> %s rejected. Acc. ratio: %2.2f"%(rej, ar)
            
            
            # Store accepted values
            for c, t,pr, acc in zip(range(self.nchains), theta, prop, accepted): #Iterates over the results of each chain
                #if not accepted repeat last value
                if not acc:
                    #Add something to the seqhist so that they all have the same length
                    if self.seqhist[c] == []:
                        self.seqhist[c].append(t)
                    else:
                        self.seqhist[c].append(self.seqhist[c][-1])
                else:
                    self.history[j, :] = t
                    self.seqhist[c].append(t)
                    try:
                        self.phi[j] = pr[0] if self.t==1 else [tuple(point) for point in pr]
                        ptheta[j] = tuple(t)
                    except IndexError:
                        print "index error",  j,  self.phi.shape
                    self.liklist.append(liks[c])
                    if j == self.samples+self.burnin:break
                    j += 1 #update accepted samples counter
            
            # Remove Outlier Chains
            if j>0 and j < self.burnin:
                outl = self._det_outlier_chains(j)
                imax = pps.tolist().index(pps.max())
                for n, c in enumerate(outl):
                    if c:
                        theta[n] = theta[imax]
                        prop[n] = prop[imax]
                        pps [n] = pps[imax]
                        liks[n] = liks[imax]
            
            el =time.time()-t0
            if int(el)%10 ==0 and el>1 and j>100:#j%100 == 0 and j>0:
                if self.trace_acceptance:
                    print "++>Acc. %s out of %s. Acc. ratio: %1.3f"%(j,i, ar)
                    self._watch_chain(j)
                if self.trace_convergence:
                    print "++> Likvar: %s\nBest run Likelihood:%s"%(self.likvariance, np.max(self.liklist) )
                t0 = time.time()
#            print "%s\r"%j
            last_pps = pps
            #last_liks = last_liks
            last_prop = prop
            last_theta = theta
            ar = (i-rej)/float(i)
            if self.meld.verbose ==2 and j > 10:
                #print len(self.liklist),j
                self.meld.current_plot(self.phi, self.data, self.best_prop_index, step=j)
        self.term_pool()
        self.meld.post_theta = ptheta[self.burnin:]
        self.meld.post_phi = self.phi[self.burnin:]
        self.meld.post_theta = ptheta#self._imp_sample(self.meld.L,ptheta,liklist)
        self.meld.likmax = max(self.liklist)
        self.meld.DIC = self.DIC
        print "Total steps(i): ",i,"rej:",rej, "j:",j
        print ">>> Acceptance rate: %1.3f"%ar
        self.shut_down('Finished normally.')
        return 1

if __name__ == "__main__":
    pass
