"""
Parameter estimation and series forcasting based on simulated data with moving window.
Stochastic model
"""
#
# Copyright 2009- by Flávio Codeço Coelho
# License gpl v3
#

from __future__ import division
from BIP.Bayes.Melding import Meld
from BIP.SDE.gillespie import Model
from scipy.integrate import odeint
import pylab as P
from numpy import array, arange,log, sqrt, ceil, floor, mean, median, diff, concatenate,cumsum,zeros,nan_to_num
from scipy import stats
from itertools import cycle
import copy
import time
import cPickle,glob
from pylab import csv2rec, rec2csv, movavg
from plotmeld import *
from BIP.Viz.realtime import RTplot
#from stoch_seir import runModel


mu = 0.0 #birth and death rate.FIXED
beta = 0.00058 #Transmission rate
eta = .5 #infectivity of asymptomatic infections relative to clinical ones. FIXED
epsilon = .1 #latency period 
alpha = .2 #Probability of developing clinical influenza symptoms
sigma = .5 #reduced risk of re-infection after recovery
tau = .01 #infectious period. FIXED
# Initial conditions
global inits,tf
tf= 140
inits = array([490,0,10,0,0])
pars = [beta,alpha,sigma]


# propensity functions
def f1(r,inits):return r[0]*inits[0]*(inits[2]+inits[3])#S->E
def f2(r,inits):return r[1]*inits[1]#E->I
def f3(r,inits):return r[3]*inits[2]#I->R
def f4(r,inits):return r[2]*inits[1]#E->A
def f5(r,inits):return r[4]*inits[3]#A->R

def runModel(*theta):
    global tf,inits
    step = 1
    #setting parameters
    beta,alpha,sigma = theta[:3]
    vnames = ['S','E','I','A','R']
    #rates: b,ki,ka,ri,ra
    #r = (0.001, 0.1, 0.1, 0.01, 0.01)
    r = (beta, (alpha)*epsilon, (1-alpha)*epsilon, tau, tau)
    #print r,inits
    # propensity functions
    propf = (f1,f2,f3,f4,f5)

    tmat = array([[-1, 0, 0, 0, 0],
                  [ 1,-1, 0,-1, 0],
                  [ 0, 1,-1, 0, 0],
                  [ 0, 0, 0, 1,-1],
                  [ 0, 0, 1, 0, 1]
                ])
    M=Model(vnames=vnames,rates = r,inits=inits,tmat=tmat,propensity=propf)
    #t0 = time.time()
    M.run(tmax=tf,reps=50,viz=False,serial=True)
    t,series,steps = M.getStats()
    ser = series.mean(axis=0)
    #print series.shape
    return ser

    
def doInference(K,L,dur,prior={},data=None):
    global tf
    Me = Meld(K,L,model=runModel, ntheta=3,nphi=5, verbose=False,viz=True)
    if prior['theta'] and prior['phi']:
        Me.setThetaFromData(['beta','alpha','sigma'],prior['theta'],[(0,1),(.001,1),(0,1)])
        Me.setPhiFromData(['S','E','I','A','R'],prior['phi'],[(0,500)]*5)
    else:
        print "===> First"
        Me.setTheta(['beta','alpha','sigma'],[stats.uniform]*3,[(0.00001,.0006),(.01,.5),(0,.02)])
        Me.setPhi(['S','E','I','A','R'],[stats.uniform]*5,[(0,500),(0,500),(0,500),(0,500),(0,500)],[(0,500)]*5)
    data={'S':data[:,0],'E':data[:,1],'I':data[:,2],'A':data[:,3],'R':data[:,4]}
    succ=0
    att = 1
    while not succ: #run sir Until is able to get a fit
        print 'attempt #',att
        succ = Me.sir(data=data,tau=100000,nopool=True,t=tf)
        att += 1
    pt,series = Me.getPosteriors(t=tf)
    pp = series[:,-1]
    # Getting predictions for following week
    # Set initial values according to data.
    # Find the simulated prevalence that best matches real value at the end of the week.
    diff = abs(pp.I-(data['I'][-1])) 
    initind = diff.tolist().index(min(diff))
    inits = list(pp[initind])
    #Adjust initial prevalence(I) by exchanging mass with the recovered compartment
    inits[2] -= pp.I[initind]-data['I'][-1] #adjusting incidence
    inits[-2] += pp.I[initind]-data['I'][-1]#adjusting Recovered
    tf = dur
    pt,predseries, = Me.getPosteriors(t=tf)
    return pt,pp,series,predseries



def testWSimData(pars):
    K = 2000
    L = 1000
    global inits,tf
    hst = RTplot()
    hsp = RTplot()
    ser = RTplot()
    weeks = 20 #number of weeks to run
    tf = 7*weeks
    subnot = stats.uniform(.8,.1).rvs((tf,5))
    beta,alpha,sigma =pars
    duration = 7 #size of the fitting window: a week
    y = runModel(*[beta,alpha,sigma])
    ser.plotlines(y.T,['S','E','I','A','R'])
    #y *=subnot #apply noise due to sub-notification
    start = time.time()
    
    prior = {'theta':[],'phi':[]}
    for w in range(weeks):
        print '==> Week # %s!'%w
        tf=7 #time frame (7 days) a week long
        if w > 0:
            inits = y[w*tf]
        #print y.shape
        #print y[w*7:w*7+7]
        pt,pp,series,predseries=doInference(K,L,duration,prior=prior,data=y[w*tf:w*tf+tf])
        f = open('week_%s'%w,'w')
        d = {'date':range(1,y.shape[0]+1),'S':y[:,0],'E':y[:,1],'I':y[:,2],'A':y[:,3],'R':y[:,4]}
        cPickle.dump((pt,pp,series, d,predseries),f)
        f.close()
        print "saved!"
        prior = {'theta':[],'phi':[]}
        for n in pt.dtype.names:
            prior['theta'].append(pt[n])
        
        for n in pp.dtype.names:
            prior['phi'].append(pp[n])
        hst.plothist(data=array(prior['theta']),names=list(pt.dtype.names),title='Week %s'%w)
        hsp.plothist(data=array(prior['phi']),names=list(pp.dtype.names),title='Week %s'%w)
        hst.clearFig()
        hsp.clearFig()
    print "time: %s seconds"%(time.time()-start)

if __name__== "__main__":
    testWSimData(pars)
    plotWeekly(data=False)
