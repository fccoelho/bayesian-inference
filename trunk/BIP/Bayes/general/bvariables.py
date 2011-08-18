"""
This module implements classes to represent an arbitrary Bayesian random variable.

This is experimental code! do not use for serious applications!

"""
# copyright 2007 Flavio Codeco Coelho
# Licensed under GPL v3
from BIP.Bayes import like
from numpy import arange
from numpy import array
from numpy import compress
from numpy import exp
from numpy import greater
from numpy import less
from numpy import ones
from numpy import searchsorted
from numpy import sqrt
import pylab as P
from scipy import stats
import sys

from BIP.Viz.asciihist import Histogram

__docformat__ = "restructuredtext en"
## Conjugate prior list: distribution types which have supported conjugate prior 
discrete_conjugate_priors = {
    'Bernoulli':'Beta',
    'Binomial':'Beta',
    'Negative Binomial':'Beta',
    'Poisson':'Gamma',
    'Multinomial':'Dirichlet',
    'Geometric':'Beta'
    }
continuous_conjugate_priors = {
    'Uniform':'Pareto',
    'Exponential':'Gamma',
    'Normal':'Normal',
    'Pareto':'Gamma',
    'Gamma':'Gamma',
    'Inverse Gamma':'Gamma'
    }


def BayesVar(priortype, pars, range, resolution=1024):
    """
    Factory function for continuous and discrete variables
    """
    if isinstance(priortype, stats.rv_continuous):
        return __BayesC(priortype, pars, range, resolution)
    if isinstance(disttype, stats.rv_discrete):
        return __BayesD(priortype, pars, range, resolution)


class _BayesVar(object):
    """
    Bayesian random variate.
    """
    def __init__(self, disttype, pars, rang, resolution=1024):
        '''
        Initializes random variable.

        :parameters:
            - `disttype`: must be a valid RNG class from scipy.stats
            - `pars`: are the parameters of the distribution.
            - `rang`: range of the variable support.
            - `resolution`: resolution of the support.
        '''

        self.distn = disttype.name
        self._flavorize(disttype(*pars), disttype)
        self.pars = pars
        self.rang = rang
        self.res = (rang[1]-rang[0]) * 1. / resolution
        self.likelihood = None
        self.data = None
        self.posterior = array([])

    def __str__(self):
        '''
        :Return:
            ascii histogram of the variable
        '''
        if self.posterior.any():
            d = self.posterior
        else:
            d = self.get_posterior_sample(200000)
        name = self.distn + self.pars.__str__()
        h = Histogram(d, bins=10)
        return name + '\n' + h.vertical()

    def _flavorize(self, pt, ptbase):
        '''
        Add methods from distribution type
        '''
        self.cdf = pt.cdf
        self.isf = pt.isf
        if isinstance(ptbase, stats.rv_continuous):
            self.pdf = pt.pdf
        elif isinstance(ptbase, stats.rv_discrete):
            self.pdf = pt.pmf
        else: sys.exit('Invalid distribution object')
        self.ppf = pt.ppf
        self.rvs = pt.rvs
    def _update(self, model):
        """
        Calculate likelihood function
        """
        if self.data != None:
            d = self.data
            sc = self.pars[1]
            m = self.rang[0]
            M = self.rang[1]
            step = self.res
            likefun = self._likelihood(model)#returns log-likelihood function
            lik = exp(array([likefun((d, i, d.var())) for i in arange(m, M, step)]))
            self.likelihood = lik / sum(lik)

    def add_data(self, data, model):
        """
        Updates variable with information from dataset

        :Parameters:
            - `data`: sequence of numbers
            - `model`: probabilistic model underlying data
        """
        self.data = array(data) 
        self._update(model.dist.name)

    def get_prior_sample(self, n):
        '''
        Returns a sample from the prior distribution

        :Parameters:
            - `n`: Sample size.
        '''
        return self.rvs(size=n)
    
    def get_prior_dist(self):
        """
        Returns the prior PDF.
        """
        return self.pdf(arange(self.rang[0], self.rang[1], self.res))
        
    def get_posterior_sample(self, n):
        """
        Return a sample of the posterior distribution.
        Uses SIR algorithm.

        :Parameters:
            - `n`: Sample size.
        """
        if self.posterior.any():# Use last posterior as prior
            k = stats.kde.gaussian_kde(self.posterior)
            s = k.resample(n)
        else:
            s = self.get_prior_sample(n)
        if self.data != None:
            m = self.rang[0]
            M = self.rang[1]
            step = self.res
            supp = arange(m, M, step)#support
            s = compress(less(s.ravel(), M) & greater(s.ravel(), m), s)#removing out-of-range samples
            d = stats.uniform.rvs(loc=0, scale=1, size=len(s))#Uniform 0-1 samples
            w = self.pdf(supp) * self.likelihood
            w = w / sum(w) #normalizing weights
            sx = searchsorted(supp, s)
            w = w[sx-1]#search sorted returns 1-based binlist
            post = compress(d < w, s)
            self.posterior = post
            return post
        else:
            return array([])

    def _likelihood(self, dname):
        '''
        Defines parametric family  of the likelihood function.
        Returns likelihood function.

        :Parameters:
            - `dname`: must be a string.
        :Return:
            lambda function to calculate  the likelihood.
        '''
        like_funs = {
            'norm': lambda(x):like.Normal(x[0], x[1], 1. / x[2]),
            'expon': lambda(x):(1. / x[2]) ** x[0].size * exp(-(1. / x[2]) * sum(x[0])),
            'beta': lambda(x):like.Beta(x[0], x[1], x[2]),
            'uniform': lambda(x): like.Uniform(x[0], x[1]-2 * sqrt(x[2]), x[1] + 2 * sqrt(x[2]))
        }
        return like_funs[dname]
        #TODO: expand for more distribution types

        
    def _post_from_conjugate(self, dname, * pars):
        '''
        Returns posterior distribution function using conjugate prior theory
        '''
        if not self.data:
            return
        if dname == 'bernoulli':
            pdist = stats.beta(pars[0])
            # TODO: finish this

class __BayesC(_BayesVar, stats.rv_continuous):
    def __init__(self, priortype, pars, range, resolution=512):
        _BayesVar.__init__(self, priortype, pars, range, resolution)

class __BayesD(_BayesVar, stats.rv_discrete):
    def __init__(self, priortype, pars, range, resolution=512):
        _BayesVar.__init__(self, priortype, pars, range, resolution)

if __name__ == "__main__":
    #bv = BayesVar(stats.norm,(3,1),range=(0,5))
    bv = BayesVar(stats.norm, (3, 1), range=(0, 5), resolution=1000)
    data = stats.uniform(1, 3).rvs(500)
    bv.add_data(data, stats.uniform(1, 3))
    print bv
    p = bv.get_posterior_sample(200000)
    print bv
    P.plot(arange(bv.rang[0], bv.rang[1], bv.res), bv.likelihood / max(bv.likelihood), 'ro', lw=2)
    P.plot(arange(bv.rang[0], bv.rang[1], bv.res), bv.get_prior_dist(), 'g+', lw=2)
    print p
    P.hist(p, normed=1)
    P.legend(['Likelihood', 'Prior', 'Posterior'])
    P.title('Bayesian inference')
    P.savefig('bayesvar.png', dpi=400)
    P.show()
