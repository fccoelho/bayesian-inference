#-*- coding:utf-8 -*-
"""
Parameter estimation and series forcasting based on simulated data with moving window.
Stochastic model
"""
#
# Copyright 2009- by Flávio Codeço Coelho
# License gpl v3
#
from BIP.SDE.gillespie import Model
from BIP.Bayes.Melding import FitModel
import numpy as np
from scipy import stats as st

mu = 0.0  #birth and death rate.FIXED
beta = 0.00058  #Transmission rate
eta = .5  #infectivity of asymptomatic infections relative to clinical ones. FIXED
epsilon = .1  #latency period
alpha = .2  #Probability of developing clinical influenza symptoms
sigma = .5  #reduced risk of re-infection after recovery
tau = .01  #infectious period. FIXED
# Initial conditions
global inits, tf
tf = 140
inits = [490, 0, 10, 0, 0]
pars = [beta, alpha, sigma]


# propensity functions
def f1(r, inits): return r[0] * inits[0] * (inits[2] + inits[3])  #S->E


def f2(r, inits): return r[1] * inits[1]  #E->I


def f3(r, inits): return r[3] * inits[2]  #I->R


def f4(r, inits): return r[2] * inits[1]  #E->A


def f5(r, inits): return r[4] * inits[3]  #A->R


def runModel(theta):
    global tf, inits
    step = 1
    #setting parameters
    beta, alpha, sigma = theta[:3]
    vnames = ['S', 'E', 'I', 'A', 'R']
    #rates: b,ki,ka,ri,ra
    #r = (0.001, 0.1, 0.1, 0.01, 0.01)
    r = (beta, alpha * epsilon, (1 - alpha) * epsilon, tau, tau)
    #print r,inits
    # propensity functions
    propf = (f1, f2, f3, f4, f5)

    tmat = np.array([[-1, 0, 0, 0, 0],
                     [1, -1, 0, -1, 0],
                     [0, 1, -1, 0, 0],
                     [0, 0, 0, 1, -1],
                     [0, 0, 1, 0, 1]
    ])
    M = Model(vnames=vnames, rates=r, inits=inits, tmat=tmat, propensity=propf)
    #t0 = time.time()
    M.run(tmax=tf, reps=1, viz=0, serial=True)
    t, series, steps, events = M.getStats()
    ser = st.nanmean(series, axis=0)
    #print series.shape
    return ser


d = runModel([beta, alpha, sigma])

dt = {'S': d[:, 0], 'E': d[:, 1], 'I': d[:, 2], 'A': d[:, 3], 'R': d[:, 4]}
F = FitModel(300, runModel, inits, tf, ['beta', 'alpha', 'sigma'], ['S', 'E', 'I', 'A', 'R'],
             wl=7, nw=20, verbose=True, burnin=100)
F.set_priors(tdists=[st.uniform] * 3, tpars=[(0.00001, .0006), (.1, .5), (0.0006, 1)],
             tlims=[(0, .001), (.001, 1), (0, 1)],
             pdists=[st.uniform] * 5, ppars=[(0, 500)] * 5, plims=[(0, 500)] * 5)

F.run(dt, 'MCMC', likvar=1e1, pool=0, monitor=[])
#~ print F.optimize(data=dt,p0=[0.1,.5,.1], optimizer='oo',tol=1e-55, verbose=1, plot=1)
#==Uncomment the line below to see plots of the results
F.plot_results()
