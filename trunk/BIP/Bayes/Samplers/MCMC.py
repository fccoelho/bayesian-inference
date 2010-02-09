# To change this template, choose Tools | Templates
# and open the template in the editor.
"""
Module implementing MCMC samplers 

    - Metropolis: Standard random Walk Metropolis Hastings sampler
    - Dream: DiffeRential Evolution Adaptive Markov chain sampler
"""

__author__="fccoelho"
__date__ ="$09/12/2009 10:44:11$"
__docformat__ = "restructuredtext en"

import numpy as np
from multiprocessing import Pool,  Process
from multiprocessing.managers import BaseManager
import scipy.stats as st
import sys
import xmlrpclib
from BIP.Viz.realtime import RTplot


class _Sampler(object):
    '''
    Base classe for all samplers
    Holds common logic and 
    '''
    _po = None
    _dimensions = None #cache for dimensions
    trace_acceptance = False
    trace_convergence = False
    _R = np.inf #Gelman Rubin Convergence
    def __init__(self,  parpriors=[],  parnames = []):
        self.parpriors = parpriors
        self.parnames = parnames

    @property
    def dimensions(self):
        if not self._dimensions:
            self._dimensions = len(self.parpriors)
        return self._dimensions
    @property
    def po(self):
        '''
        Pool of processes for parallel execution of tasks
        Remember to call self.po.close() and self.po.join() when done.
        '''
        if self._po == None:
            self._po = Pool()
        else:
            if self_po._state:
                self._po = Pool() #Pool has been terminated
        return self._po

    def term_pool(self):
        if self._po == None:
            return
        if not self_po._state: #Pool needs terminating
            self_po.close()
            self._po.join()
            self._po = None
    def gr_R(self):
        return self._R
    def gr_convergence(self, sequences, relevantHistoryEnd, relevantHistoryStart):
        """
        """
        pass

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
        self.phi = np.recarray((self.samples+self.burnin,t),formats=['f8']*self.meld.nphi, names = self.meld.phi.dtype.names)
        self.scaling_factor = (2.4**2)/self.dimensions
        self.e = 1e-20
        self.history = []
        if kwargs:
            for k, v in kwargs.iteritems():
                exec('self.%s = %s'%(k, v))
        
        self.chp = RTplot()

    def _propose(self, step=0):
        """
        Generates proposals
        """
        if step <= 1: 
            #sample from the priors
            while 1:
                theta = [self.parpriors[dist]() for dist in self.parnames]
                if sum ([int(np.greater(t, self.parlimits[i][0]) and np.less(t, self.parlimits[i][1])) for i, t in enumerate(theta)]) == self.dimensions:
                    break
            self.lastcv = None
        else:
            #use gaussian proposal
            if (self.lastcv==None) or (step%3==0): #recalculate covariance matrix only every ten steps
                cv = self.scaling_factor*st.cov(np.array(self.history))+self.scaling_factor*self.e*np.identity(self.dimensions)
                self.lastcv = cv
            else:
                cv = self.lastcv
            while 1:
                theta = np.random.multivariate_normal(self.history[-1],cv, size=1).tolist()[0]
                if sum ([int(np.greater(t, self.parlimits[i][0]) and np.less(t, self.parlimits[i][1])) for i, t in enumerate(theta)]) == self.dimensions:
                    break
        prop = self.meld.model_as_ra(theta)[:self.t]
#        if self.t == 1:
#            prop=[(v,) for v in prop]

        return theta,prop

    def step(self):
        """
        Does the actual sampling loop.
        """
        ptheta = np.recarray(self.samples+self.burnin,formats=['f8']*self.dimensions, names = self.parnames)
        i=0;j=0;rej=0;ar=0 #total samples,accepted samples, rejected proposals, acceptance rate
        last_lik = None
        liklist = []
        while j < self.samples+self.burnin:
            #generate proposals
            theta,prop = self._propose(step=j)
            #calculate likelihoods
            lik = self.meld._output_loglike(prop, self.data, self.likfun, self.likvariance)
            accepted = self._accept(last_lik, lik) if last_lik else 0
            # ABC fit
#            lik = self._rms_fit(prop, self.data)
#            accepted = self._accept(last_lik, lik) if last_lik else 0
            #Decide whether to accept proposal
            if last_lik == None: #on first sample
                last_lik = lik
                continue
            if not accepted:
                rej +=1 #adjust rejection counter
                i +=1
                if rej%100 == 0: 
                    #self._tune_likvar(ar)
                    print "-->%s Rejected. Last Proposal: %s"%(rej,theta)
                    print j," Accepted. likvar: %s"%self.likvariance
                continue
            self.history.append(theta)
            self.phi[j] = prop[0] if self.t==1 else [tuple(point) for point in prop]
            ptheta[j] = tuple(theta)
            liklist.append(lik)
            #self._tune_likvar(ar)
            if j%100==0:
                self._watch_chain()
                if self.trace_acceptance:print "++>%s: Acc. ratio: %s"%(j, ar)
                if self.trace_convergence: print "%s: Mean Proposal: %s; STD: %s"%(j, np.array(self.history[-100:]).mean(axis=0),np.array(self.history[-100:]).std(axis=0) )
            last_lik = lik
            j += 1 #update good sample counter 
            i+=1
            ar = j/(float(j)+rej)
        self.meld.post_theta = ptheta[self.burnin:]
        self.phi = self.phi[self.burnin:]
        self.meld.post_theta = self._imp_sample(self.meld.L,ptheta,liklist)
        print "Total steps(i): ",i,"rej:",rej, "j:",j
        print ">>> Acceptance rate: %s"%ar
        self.term_pool()
        return 1
    def _tune_likvar(self, ar):
        if ar<=.1:
            self.likvariance *= 5
        elif ar <= .15:
            self.likvariance *= 2
        elif ar <= .2:
            self.likvariance *= 1.2
        elif ar >.7:
            self.likvariance *= 0.8
        elif ar > .8:
            self.likvariance *= 0.5
        elif ar > .9:
            self.likvariance *= 0.2
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
        
    def _accept(self,  last_lik,  lik):
        """
        Decides whether to accept a proposal
        """
        # liks are logliks
        if lik == -np.inf:#0:
            return 0
        if last_lik >-np.inf:#0:
            alpha = min( np.exp(lik-last_lik), 1)
            #alpha = min(lik-last_lik, 1)
        elif last_lik == -np.int:#0:
            alpha = 1
        else:
            return 0
            raise ValueError("Negative likelihood!?!")
#        print "last_lik, lik, alpha: ",  last_lik, lik, alpha
        if np.random.random() < alpha:
            return 1
        else:
            return 0

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

    def _watch_chain(self):
        if len(self.history)<100:
            return
        s = xmlrpclib.ServerProxy('http://localhost:9876')
        s.clearFig()
        data = np.array(self.history[-100:]).T.tolist()
        s.plotlines(data, self.parnames, "Chain Progress") 

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

class Dream(_Sampler):
    '''
    DiffeRential Evolution Adaptive Markov chain sampler
    '''
    def __init(self, samples = 1000, sampmax = 20000 , parpriors=[], parnames=[], burnIn = 100, thin = 5, convergenceCriteria = 1.1,  nCR = 3, DEpairs = 1, adaptationRate = .65, eps = 5e-6, mConvergence = False, mAccept = False):
        self.samples = samples
        self.sampmax = sampmax
        self.parpriors = parpriors
        self.parnames = parnames
        self.nChains = len(parpriors)
        self.maxChainDraws = floor(samples/self.nChains)
        self.nCR = nCR
        self.DEpairs = DEpairs
        #initialize the history arrays   
        self.combinedHistory = zeros((self.nChains * self.maxChainDraws , self.dimensions))
        self.sequenceHistories = zeros((self.nChains, self.dimensions, self.maxChainDraws))
        # initialize the temporary storage vectors
        self.currentVectors = zeros((nChains, dimensions))
        self.currentLiks = zeros(nChains)
        self.scaling_factor = 2.38/np.sqrt(2*Depairs*self.dimensions)

    def _propose(self,  step=0):
        """
        Generates proposals

        :Returns:
            - theta: list of parameters vectors of length `self.dimensions`
            - prop: list of output arrays of length `self.dimensions`
        """
        #Draw a Theta for each chain
        for i in range(self.nChains):
            if step == 0:
                theta = [self.parpriors[dist](self.dimensions) for dist in self.parnames]
                self.lastcv = None
            else:
                #draw from gaussian proposal dist
                if (self.lastcv==None) or (step%1==0): #recalculate covariance matrix only every ten steps
                    cv = self.scaling_factor*st.cov(np.array(self.history))+self.scaling_factor*self.e*np.identity(self.dimensions)
                    self.lastcv = cv
                else:
                    cv = self.lastcv
                while 1:
                    theta = np.random.multivariate_normal(self.history[-1],cv, size=1).tolist()[0]
                    if sum ([int(np.greater(t, self.parlimits[i][0]) and np.less(t, self.parlimits[i][1])) for i, t in enumerate(theta)]) == self.dimensions:
                        break
        prop = self.meld.model_as_ra(theta)[:self.t]
       

    def _get_likelihoods(prop, po = None):
        '''
        Calculate the likelihoods for each chain
        '''
        if po:
            listoliks = [po.apply_async(self.meld._output_loglike, (p, self.data, self.likfun, self.likvariance)) for p in prop]
            listoliks = [l.get() for l in listoliks]
            self.term_pool()
        else:
            listoliks = [self.meld._output_loglike(p, self.data, self.likfun, self.likvariance) for p in prop]
        return listoliks
        
    def step(self):
        """
        """
        
if __name__ == "__main__":
    pass
