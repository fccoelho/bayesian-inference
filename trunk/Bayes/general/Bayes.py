"""
This module implements classes to represent an arbitrary Bayesian random variable.

"""
# copyright 2007 Flavio Codeco Coelho
# Licensed under GPL v3
from numpy import arange,compress, array, exp
import like, sys
import pylab as P
from scipy import stats

## Conjugate prior list: distribution types which have supported conjugate prior 

conjlist = [
    'bernoulli',
    'binom',
    'nbinom', #negative binomial
    'poisson',
    'geom', #geometric
    ]

## Factory functions for continuous and discrete variables 
def Continuous(priortype,pars, range,resolution=512):
    return __BayesC(priortype,pars, range,resolution)

def Discrete(priortype,pars, range,resolution=512):
    return __BayesD(priortype,pars, range,resolution)

class _BayesVar(object):
    """
    Bayesian random variate.
    """
    def __init__(self, disttype,pars, rang,resolution=512):
        '''
        Initializes random variable.
        
        :Parameters:
        - `disttype`: must be a valid RNG from scipy.stats
        - `pars`: are the parameters of the distribution.
        '''
        self.distn = dist_type.name
        self._flavorize(disttype(*pars), disttype)
        self.pars = pars
        self.rang = rang
        self.res = (rang[1]-rang[0])*1./resolution
        self.likefun = self._Likelihood(self.distn)
        self.likelihood = None
        self.data = []
        self.posterior=array([])

    def _flavorize(self,pt, ptbase):
        '''
        Add methods from distribution type
        '''
        self.cdf = pt.cdf
        self.isf = pt.isf
        if isinstance(ptbase,stats.rv_continuous):
            self.pdf = pt.pdf
        elif isinstance(ptbase,stats.rv_discrete):
            self.pdf = pt.pmf
        else: sys.exit('Invalid distribution object')
        self.ppf = pt.ppf
        self.rvs = pt.rvs
    def _update(self):
        """
        Calculate likelihood function
        """
        if self.data:
            d = self.data[-1]
            sc = self.pars[1]
            m = self.rang[0]
            M = self.rang[1]
            step = self.res
            #self.likefun returns log-likelihood
            lik = exp(array([self.likefun((d,i,sc)) for i in arange(m,M,step)]))
            self.likelihood = lik/sum(lik)

    def addData(self, data = []):
        """
        Adds dataset to variable's data store
        """
        self.data.append(array(data))
        self._update()

    def getPriorSample(self,n):
        '''
        Returns a sample from the prior distribution
        '''
        return self.rvs(size=n)
    
    def getPriorDist(self):
        """
        Returns the prior PDF.
        """
        return self.pdf(arange(self.rang[0],self.rang[1],self.res))
        
    def getPosteriorSample(self, n):
        """
        Return a sample of the posterior distribution.
        Uses SIR algorithm.
        """
        if self.posterior.any():# Use last posterior as prior
            k= stats.kde.gausian_kde(self.posterior)
            s= k.resample(n)
        else:
            s = self.getPriorSample(n)
        if self.data:
            m = self.rang[0]
            M = self.rang[1]
            step = self.res
            supp = arange(m,M,step)#support
            s = compress(less(s.ravel(),M) & greater(s.ravel(),m),s)#removing out-of-range samples
            d = stats.uniform.rvs(loc=0,scale=1,size=len(s))#Uniform 0-1 samples
            w = self.pdf(supp)*self.likelihood
            w = w/sum(w) #normalizing weights
            sx = searchsorted(supp,s)
            w = w[sx-1]#search sorted returns 1-based binlist
            post = compress(d<w,s)
            self.posterior = post
            return post
        else:
            return array([])

    def _Likelihood(self,typ):
        '''
        Defines parametric family  of the likelihood function.
        Returns likelihood function.
        typ must be a string.
        '''
        if typ == 'norm':
            return lambda(x):like.Normal(x[0],x[1],1./x[2])
        elif typ == 'expon':
            return lambda(x):(1./x[2])**x[0].size*exp(-(1./x[2])*sum(x[0]))
        elif typ == 'beta':
            return lambda(x):like.Beta(x[0],x[1],x[2])
        
    def _postFromConjugate(dname,*pars):
        '''
        Returns posterior distribution function using conjugate prior theory
        '''
        if not self.data:
            return
        if dname == 'bernoulli':
            pdist = stats.beta(pars[0])
            # TODO: finish this

class __BayesC(_BayesVar, stats.rv_continuous):
    def __init__(self, priortype,pars, range,resolution=512):
        _BayesVar.__init__(self, priortype,pars, range,resolution)

class __BayesD(_BayesVar, stats.rv_discrete):
    def __init__(self, priortype,pars, range,resolution=512):
        _BayesVar.__init__(self, priortype,pars, range,resolution)
if __name__=="__main__":
    #bv = BayesVar(stats.norm,(3,1),range=(0,5))
    bv = Continuous(stats.norm,(3,1),range=(0,5), resolution=1000)
    data = ones(10)
    bv.addData(data)
    p = bv.getPosteriorSample(200000)
    P.plot(arange(bv.rang[0],bv.rang[1], bv.res),bv.likelihood/max(bv.likelihood), 'ro', lw=2)
    P.plot(arange(bv.rang[0],bv.rang[1], bv.res),bv.getPriorDist(),'g+',lw=2)
    P.hist(p, normed=1)
    P.legend(['Likelihood','Prior', 'Posterior'])
    P.title('Bayesian inference')
    P.savefig('bayesvar.png',dpi=400)
    P.show()
