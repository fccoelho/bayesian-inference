#-*- coding: utf-8 -*-
"""
Parameter estimation and series forcasting based on simulated data with moving window.
Deterministic model
"""
#
# Copyright 2009- by Flávio Codeço Coelho
# License gpl v3
#
from BIP.Bayes.Melding import FitModel
from scipy.integrate import odeint
import scipy.stats as st
import numpy as np

beta = 1 #Transmission coefficient
tau = .2 #infectious period. FIXED
tf = 36
y0 = [.999,0.001,0.0]
def model(theta):
    beta = theta[0]
    def sir(y,t):
        '''ODE model'''
        S,I,R = y
        return  [-beta*I*S, #dS/dt
                beta*I*S - tau*I, #dI/dt
                tau*I] #dR/dt
    y = odeint(sir,inits,np.arange(0,tf,1))
    return y
    
F = FitModel(500, model,y0,tf,['beta'],['S','I','R'],
            wl=36,nw=1,verbose=0,burnin=100)
F.set_priors(tdists=[st.norm],tpars=[(1.1,.2)],tlims=[(0.5,1.5)],
    pdists=[st.uniform]*3,ppars=[(0,.1),(0,.1),(.8,.2)],plims=[(0,1)]*3)
d = model([1.0]) #simulate some data
noise = st.norm(0,0.01).rvs(36)
dt = {'I':d[:,1]+noise} # add noise
F.run(dt,'MCMC',likvar=1e-5,pool=True,monitor=[])
#==Uncomment the line below to see plots of the results
F.plot_results()
