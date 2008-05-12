#-*-encoding:latin-1-*-
# Disponivel no pacote de programas como: Bayes.py
# copyright 2007 Flavio Codeco Coelho
# Licensed under GPL v3
from numpy import * 
import like, sys
import pylab as P
from scipy.stats import *    

class BayesVar(object):
    """
    Bayesian random variate.
    """
    def __init__(self, priortype,pars, range,resolution=512):
        '''
        Inicializa variável aleatória.
        Adquire métodos da classe  priortype.
        priortype deve ser um RNG de scipy.stats
        pars são os parâmetros da priori.
        '''
        self.priorn = priortype.name
        self._flavorize(priortype(*pars), priortype)
        self.pars = pars
        self.range = range
        self.res = (range[1]-range[0])*1./resolution
        self.likefun = self._Likelihood(self.priorn)
        self.likelihood = None
        self.data = []
        self.posterior=array([])

    def _flavorize(self,pt, ptbase):
        '''
        add methods from distribution type
        '''
        self.cdf = pt.cdf
        self.isf = pt.isf
        if isinstance(ptbase,rv_continuous):
            self.pdf = pt.pdf
        elif isinstance(ptbase,rv_discrete):
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
            m = self.range[0]
            M = self.range[1]
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
        return self.pdf(arange(self.range[0],self.range[1],self.res))
    def getPosteriorSample(self, n):
        """
        Return a sample of the posterior distribution.
        """
        if self.posterior.any():# Use last posterior as prior
            k= kde.gausian_kde(self.posterior)
            s= k.resample(n)
        else:
            s = self.getPriorSample(n)
        if self.data:
            m = self.range[0]
            M = self.range[1]
            step = self.res
            supp = arange(m,M,step)#support
            s = compress(less(s.ravel(),M) & greater(s.ravel(),m),s)#removing out-of-range samples
            d = uniform.rvs(loc=0,scale=1,size=len(s))#Uniform 0-1 samples
            w = self.pdf(supp)*self.likelihood
            w = w/sum(w) #normalizing weights
            sx = searchsorted(supp,s)
            w = w[sx-1]#search sorted returns 1-based binlist
            post = compress(d<w,s)
            self.posterior = post
            return post
        else:
            return []

    def _Likelihood(self,typ):
        '''
        Define familia paramétrica da verossimilhança.
        Retorna função de verossimilhança.
        typ deve ser uma string.
        '''
        if typ == 'norm':
            return lambda(x):like.Normal(x[0],x[1],1./x[2])
        elif typ == 'expon':
            return lambda(x):(1./x[2])**x[0].size*exp(-(1./x[2])*sum(x[0]))
        elif typ == 'beta':
            return lambda(x):like.Beta(x[0],x[1],x[2])

if __name__=="__main__":
    bv = BayesVar(norm,(3,1),range=(0,5))
    data = ones(20)
    bv.addData(data)
    p = bv.getPosteriorSample(200000)
    P.plot(arange(bv.range[0],bv.range[1], bv.res),bv.likelihood, 'ro', lw=2)
    P.plot(arange(bv.range[0],bv.range[1], bv.res),bv.getPriorDist(),'g+',lw=2)
    P.hist(p, normed=1)
    P.legend(['Likelihood','Prior'])
    P.title('Bayesian inference')
    P.savefig('bayesvar.png',dpi=400)
    P.show()
