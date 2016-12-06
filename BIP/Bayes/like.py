# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
# Name:        like.py
# Project:  BayesianInference
# Purpose:     log-likelihood functions
#
# Author:      Flávio Codeço Coelho <fccoelho@gmail.com>
#
# Created:     2008-11-26
# Copyright:   (c) 2008 by the Author
# Licence:     GPL v3
#-----------------------------------------------------------------------------
from __future__ import absolute_import
__docformat__ = "restructuredtext en"
from six.moves import range
import scipy
import numpy as np
from scipy.misc import factorial
from scipy.special import gammaln
from numpy import array, searchsorted, log, random, pi, sum, inf
from functools import lru_cache



def Categor(x, hist):
    """
    Categorical Log-likelihood
    generalization of a Bernoulli process for variables with any constant
    number of discrete values.
    
    :Parameters:
        - `x`: data vector (list)
        - `hist`: tuple (prob,classes) classes contain the superior limit of the histogram classes
    
    >>> Categor([1],([.3,.7],[0,1]))
    -0.356674943939
    """
    like = 0.0
    x = array(x)
    prob = array(hist[0])
    sup = array(hist[1])
    ind = searchsorted(sup, x)
    like += sum(log(prob[ind]))
    return like


def Uniform(x, xmin, xmax):
    """
    Uniform Log-likelihood

    :Parameters:
        - `x`: data vector(list)
        - `min`: lower limit of the distribution
        - `max`: upper limit of the distribution

    >>> Uniform([1.1,2.3,3.4,4],0,5)
    -6.4377516497364011
    >>> Uniform([1.1,2.3,3.4,6],0,5)
    -inf
    """
    assert xmax > xmin
    x = array(x)
    like = 0.0
    p = 1. / xmax - xmin
    support = ((x > xmin) & (x <= xmax))*1.
    support[np.in1d(support, 0)] = inf
    like = support * log(p)
    # print(like)



    return np.sum(like)

@lru_cache(maxsize=1024)
def Normal(x, mu, tau):
    """
    Normal Log-like 
    
    :Parameters:
        -  `mu`: mean
        -  `tau`: precision (1/variance)
    
    >>> Normal((0,),0,1)
    -0.918938533205
    """
    x = array(x)
    n = x.size
    like = sum(-0.5 * tau * (x - mu) ** 2)
    like += n * 0.5 * log(0.5 * tau / pi)
    return like

@lru_cache(maxsize=1024)
def find_best_tau(x, mu):
    """
    returns the value of tau which maximizes normal loglik for a fixed (x,mu)
    :param x:
    :param mu:
    """

    tau = 1. / (mu + 1) if mu == 0 else 1. / mu  #starting point
    ll = Normal(x, mu, tau)
    i = 0;
    j = 0
    while i < 1000 and j < 100000:
        taun = tau + random.normal()
        l = Normal(x, mu, taun)
        if l > ll:
            tau = taun
            ll = l
            i += 1
        j += 1
    return tau


def Lognormal(x, mu, tau):
    """
    Lognormal Log-likelihood
    
    :Parameters:
        -  `mu`: mean
        -  `tau`: precision (1/sd)
    
    >>> Lognormal((0.5,1,1.2),0,0.5)
    -3.15728720569
    """
    x = array(x)
    n = x.size
    like = n * 0.5 * (log(tau) - log(2.0 * pi)) + sum(0.5 * tau * (log(x) - mu) ** 2 - log(x))
    return -like


def Poisson(x, mu):
    """
    Poisson Log-Likelihood function
    :param x: vector of data
    :param mu: mean

    >>> Poisson([2],2)
    -1.30685281944
    """
    x = array(x)
    sumx = sum(x * log(mu) - mu)
    sumfact = sum(log(factorial(x)))
    like = sumx - sumfact
    return like


def Negbin(x, r, p):
    """
    Negative Binomial Log-Likelihood 

    :param x: data
    :param r:
    :param p:

    >>> Negbin([2,3],6,0.3)
    -9.16117424315
    """
    x = array(x)
    like = sum(r * log(p) + x * log(1 - p) + log(factorial(x + r - 1)) - log(factorial(x)) - log(
        factorial(r - 1)))
    return like


def Binomial(x, n, p):
    """
    Binomial Log-Likelihood 

    :param x: data
    :param n:
    :param p:

    >>> Binomial([2,3],6,0.3)
    -2.81280615454
    """
    x = array(x)
    like = sum(x * log(p) + (n - x) * log(1. - p) + log(factorial(n)) - log(factorial(x)) - log(
        factorial(n - x)))
    return like


def Weibull(x, alpha, beta):
    """
    Log-Like Weibull

    :param x: data
    :param alpha:
    :param beta:

    >>> Weibull([2,1,0.3,.5,1.7],1.5,3)
    -7.811955373
    """
    x = array(x)
    beta = float(beta)
    n = x.size
    #Normalizing constant
    like = n * (log(alpha) - alpha * log(beta))
    # Kernel of the distribution
    like += sum((alpha - 1) * log(x) - (x / beta) ** alpha)
    return like


def Bernoulli(x, p):
    """
    Log-Like Bernoulli

    :param x: data
    :param p: probability

    >>> Bernoulli([0,1,1,1,0,0,1,1],0.5)
    -5.54517744448
    """
    x = array(x)
    like = sum(x * log(p) + (1 - x) * log(1. - p))
    return like


def Gamma(x, alpha, beta):
    """
    Log-Like Gamma

    :param x: data
    :param alpha:
    :param beta:

    >>> Gamma([2,3,7,6,4],2,2)
    -11.015748357
    """
    x = array(x)
    beta = float(beta)
    n = x.size
    #Normalizing constant
    like = -n * (gammaln(alpha) + alpha * log(beta))
    # Kernel of the distribution
    like += sum((alpha - 1.0) * log(x) - x / beta)
    return like


def Beta(x, a, b):
    """
    Log-Like Beta

    :param x: data
    :param a:
    :param b:

    >>> Beta([.2,.3,.7,.6,.4],2,5)
    -0.434845728904
    """
    x = array(x)
    n = x.size
    #Normalizing constant
    like = n * (gammaln(a + b) - gammaln(a) - gammaln(b))
    # Kernel of the distribution
    like += sum((a - 1.0) * log(x) + (b - 1.0) * log(1.0 - x))
    return like


def Simple(x, w, a, start=0):
    """
    Find out what it is. ;-)
    """
    m = len(a)
    n = len(x)
    like = 0.0
    s = sum(a * (x / w) ** (2 * list(range(n))))
    like += log(1 + s)
    return like


if __name__ == "__main__":
    import doctest

    doctest.testmod(verbose=True)
