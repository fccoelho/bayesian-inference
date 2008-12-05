#!/usr/bin/python
# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
# Name:        lhs.py
# Project:  Bayesian-Inference
# Purpose:     
#
# Author:      Flávio Codeço Coelho<fccoelho@gmail.com>
#
# Created:     2008-11-26
# Copyright:   (c) 2008 by the Author
# Licence:     GPL
#-----------------------------------------------------------------------------
__docformat__ = "restructuredtext en"
from pylab import plot, figure,hist,show, savefig, legend
import scipy.stats as stats
import numpy

def lhs(dist, parms, n=100):
    '''
    Latin Hypercube sampling of any distrbution.
    dist is is a scipy.stats random number generator 
    such as stats.norm, stats.beta, etc
    parms is a tuple with the parameters needed for 
    the specified distribution.
    '''
    perc = numpy.arange(0,1.,1./n)
    numpy.random.shuffle(perc)
    smp = [stats.uniform(i,1./n).rvs()[0] for i in perc]
    v = dist(*parms).ppf(smp)
    return v
            
if __name__=='__main__':
    dist = stats.norm
    #dist = stats.beta
    pars = (50,1)
    #pars = (1,5) #beta
    c=lhs(dist, pars,20)
    hist(c,normed=1, label='LHS sample')
    n = dist(*pars).rvs(size=20)
    hist(n.ravel(),facecolor='r',alpha =0.3,normed=1, label='Regular sample')
    plot(numpy.arange(min(min(c),min(n)),max(max(c),max(n)),.1),dist(*pars).pdf(numpy.arange(min(min(c),min(n)),max(max(c),max(n)),.1)),label='PDF')
    legend()#['LHS sample','regular sample', 'PDF'])
    #savefig('lhs.png',dpi=400)
    show()
    
#TODO: Add lhs from density, and lhs from sample
#TODO: Add correlated multiple lhs sampling