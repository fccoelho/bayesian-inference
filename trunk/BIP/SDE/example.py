# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
# Name:        example.py
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
"""
Example of an SEIR model with two Infectious classes: subclinical(Is) and clinical(Ic)
        Is
       /  \
S -> E     R
       \  /
        Ic
        
States:
S: Susceptible
E: Exposed
Is: Infectious subclinical
Ic: Infectious clinical
R: Recovered

Transition rates:
b,ks,kc,rs,rc = (0.001, 0.1, 0.1, 0.01, .01)
Transitions:
S -> E : b*S*(Is+Ic)
E -> Is : ks*E
E -> Ic : kc*E
Is -> R : rs*Is
Ic -> R : rc*Ic

"""
from gillespie import Model
import time

from numpy import array
vnames = ['S','E','Is','Ic','R']
#rates: b,ks,kc,rs,rc

r = (0.001, 0.1, 0.1, 0.01, .01)
ini = (490,0,0,10,0)
# propensity functions
def f1(r,ini):return r[0]*ini[0]*(ini[2]+ini[3])
def f2(r,ini):return r[1]*ini[1]
def f3(r,ini):return r[2]*ini[1]
def f4(r,ini):return r[3]*ini[2]
def f5(r,ini):return r[4]*ini[3]

propf = (f1,f2,f3,f4,f5)

tmat = array([[-1, 0, 0, 0, 0],#S
              [ 1,-1,-1, 0, 0],#E
              [ 0, 1, 0,-1, 0],#Is
              [ 0, 0, 1, 0,-1],#Ic
              [ 0, 0, 0, 1, 1]])#R
M=Model(vnames=vnames,rates = r,inits=ini,tmat=tmat,propensity=propf)
t0 = time.time()
M.run(tmax=80,reps=1000,viz=False,serial=False)
print 'total time: ',time.time()-t0, ' seconds.'
t,series,steps = M.getStats()
print steps,'steps'
#print series.shape
from pylab import plot , show, legend, errorbar
#print series.shape
plot(t,series.mean(axis=0),'-o')
legend(vnames,loc=0)
show()
