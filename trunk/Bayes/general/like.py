
# log-likelihood functions
# copyright 2007 by Flavio Codeco Coelho
# License: GPL
#
import scipy
from scipy.special import gammaln
from numpy import *


def Categor(x,hist):
    """
    Categorical Log-likelihood
    x: data vector (list)
    hist: tuple (prob,classes) classes contain the superior limit of the histogram classes
    >>> Categor([1],([.3,.7],[0,1]))
    -0.356674943939
    """
    like =0.0
    x = array(x)
    prob = array(hist[0])
    sup = array(hist[1])
    ind = searchsorted(sup,x)
    like += sum(log(prob[ind]))
    return like

def Normal(x,mu,tau):
    """
    Normal Log-like 
    mu: mean
    tau: precision (1/sd)
    >>> Normal([0],0,1)
    -0.918938533205
    """
    x = array(x)
    n = x.size
    like =  sum(-0.5 * tau * (x-mu)**2)
    like += n*0.5*log(0.5*tau/pi)
    return like
    
def Lognormal(x,mu,tau):
    """
    Lognormal Log-likelihood
    mu: mean
    tau: precision (1/sd)
    >>> Lognormal([0.5,1,1.2],0,0.5)
    -3.15728720569
    """
    x = array(x)
    n = x.size
    like = n * 0.5 * (log(tau) - log(2.0*pi)) + sum(0.5*tau*(log(x)-mu)**2 - log(x))
    return like

def Poisson(x,mu):
    """
    Log-Like Poisson
    >>> Poisson([2],2)
    -1.30685281944
    """
    x=array(x)
    sumx = sum(x*log(mu)-mu)
    sumfact = sum(log(scipy.factorial(x)))
    like = sumx-sumfact
    return like

def Negbin(x,r,p):
    """
    Log-Like negative Binomial
    >>> Negbin([2,3],6,0.3)
    -9.16117424315
    """
    x = array(x)
    like = sum(r*log(p)+x*log(1-p)+log(scipy.factorial(x+r-1))-log(scipy.factorial(x))-log(scipy.factorial(r-1)))
    return like

def Binomial(x,n,p):
    """
    Log-Like Binomial
    >>> Binomial([2,3],6,0.3)
    -2.81280615454
    """
    x = array(x)
    like = sum(x*log(p)+ (n-x)*log(1.-p)+log(scipy.factorial(n))-log(scipy.factorial(x))-log(scipy.factorial(n-x)))
    return like

def Weibull(x,alpha,beta):
    """
    Log-Like Weibull
    >>> Weibull([2,1,0.3,.5,1.7],1.5,3)
    -7.811955373
    """
    x=array(x)
    beta = float(beta)
    n = x.size
    #Normalizing constant
    like = n * (log(alpha) - alpha*log(beta))
    # Kernel of the distribution
    like += sum((alpha-1) * log(x) - (x/beta)**alpha)
    return like

def Bernoulli(x,p):
    """
    Log-Like Bernoulli
    >>> Bernoulli([0,1,1,1,0,0,1,1],0.5)
    -5.54517744448
    """
    x=array(x)
    like = sum(x*log(p) + (1-x)*log(1.-p))
    return like

def Gamma(x,alpha,beta):
    """
    Log-Like Gamma
    >>> Gamma([2,3,7,6,4],2,2)
    -11.015748357
    """
    x = array(x)
    beta = float(beta)
    n = x.size
    #Normalizing constant
    like = -n * (gammaln(alpha) + alpha*log(beta))
    # Kernel of the distribution
    like += sum((alpha - 1.0)*log(x) - x/beta)
    return like
    
def Beta(x,a,b):
    """
    Log-Like Beta
    >>> Beta([.2,.3,.7,.6,.4],2,5)
    -0.434845728904
    """
    x = array(x)
    n = x.size
    #Normalizing constant
    like = n * (gammaln(a+b) - gammaln(a) - gammaln(b))
    # Kernel of the distribution
    like += sum((a-1.0)*log(x) + (b-1.0)*log(1.0-x))
    return like

def Simple(x,w,a,start=0):
    """
    find out what it is.
    not implemented.
    """
    m=len(a)
    n=len(x)
    pass
    
if __name__=="__main__":
    import doctest
    doctest.testmod(verbose=True)
