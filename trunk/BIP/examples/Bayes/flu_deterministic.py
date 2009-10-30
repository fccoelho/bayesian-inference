
"""
Inference based on simulated data with moving windows.
"""

from __future__ import division
from BIP.Bayes.Melding import Meld
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


mu = 0.0 #birth and death rate.FIXED
beta = 1.4 #Transmission coefficient
eta = .5 #infectivity of asymptomatic infections relative to clinical ones. FIXED
epsilon = .2 #latency period 
alpha = .2 #Probability of developing clinical influenza symptoms
sigma = .5 #reduced risk of re-infection after recovery
tau = .5 #infectious period. FIXED
# Initial conditions
global inits,tf
tf= 140
inits = array([.999,0,.001,0,0,0])
pars = [beta,alpha,sigma]

def runModel(*theta):
    global inits,tf
    y0 = inits
    t0 = 0
    step = 1
    #setting parameters
    beta,alpha,sigma = theta[:3]
    # setting the initial number of susceptible and recovered
    #y0[-2] = theta[5]
    #y0[0] = 1-sum(y0[1:])
    #print 'hi!', y0
    def seir(y,t):
        '''ODE model'''
        S,E,I,A,R,Ia = y
        lamb = beta*(I+eta*A)
        return  [mu - (lamb+mu)*S, #dS/dt
                lamb*S - (epsilon+mu)*E, #dE/dt
                alpha*epsilon*E + alpha*sigma*lamb*R - (tau+mu)*I, #dI/dt
                (1-alpha)*epsilon*E + (1-alpha)*sigma*lamb*R - (tau+mu)*A, #dA/dt
                tau*I + tau*A - R*(sigma*lamb+mu), #dR/dt
                alpha*epsilon*E + alpha*sigma*lamb*R #Incidencia acumulada dIa/dt
                ]
    y = odeint(seir,y0,arange(t0,tf,step))
    return y
    
def doInference(K,L,dur,prior={},data=None, method='SIR'):
    global tf,inits
    Me = Meld(K,L,model=runModel, ntheta=3,nphi=6, verbose=False,viz=False)
    if prior['theta'] and prior['phi']:
        Me.setThetaFromData([r'$\beta$',r'$\alpha$',r'$\sigma$'],prior['theta'],[(0,5),(0,1),(0,1)])
        Me.setPhiFromData(['S','E','I','A','R','Ia'],prior['phi'],[(0,1)]*6)
    else:
        print "===> First"
        Me.setTheta([r'$\beta$',r'$\alpha$',r'$\sigma$'],[stats.uniform]*3,[(1.1,.7),(.01,.5),(0.4,.2)])
        Me.setPhi(['S','E','I','A','R','Ia'],[stats.uniform]*6,[(0,1),(0,1),(0,1),(0,1),(0,1),(.0,1)],[(0,1)]*6)
    data={'S':data[:,0],'E':data[:,1],'I':data[:,2],'A':data[:,3],'R':data[:,4],'Ia':data[:,5]}
    succ=0
    att = 1
    if method == 'SIR':
        while not succ: #run sir until is able to get a fit
            print 'attempt #',att
            succ = Me.sir(data=data,variance=0.0001,nopool=True,t=tf)
            att += 1
    elif method == 'ABC':
        while not succ: #run sir until is able to get a fit
            print 'attempt #',att
            succ = Me.abcRun(data=data,fitfun=dist_calc,t=tf,nopool=True)
            att += 1
    pt,series = Me.getPosteriors(t=tf)
    pp = series[:,-1]
    # Getting predictions for following week
    # Set initial values according to data.
    # Find the simulated cumulative prevalence that best matches real value at the end of the week.
    diff = abs(pp.Ia-(data['I'][-1])) 
    initind = diff.tolist().index(min(diff))
    inits = list(pp[initind])
    #Adjust initial prevalence(I) by exchanging mass with the recovered compartment
    #inits[0] = data['S'][-1]
    #inits[1] = data['E'][-1]
    inits[2] = data['I'][-1] #adjusting incidence
    inits[3] = data['A'][-1]
    #inits[3] = inits[2]*(pp.A[initind]/pp.I[initind]) #adjusting asymptomatic incidence based on proportion of asymptomatics
    inits[-1] = data['Ia'][-1]#adjusting cumulative incidence
    #inits[-2] = 1-(sum(inits[:-2])+inits[-1])#adjusting Recovered
    #inits[0] = 1-(sum(inits[1:]))#Adjusting Susceptible
    
    pt,predseries = Me.getPosteriors(t=dur)
    return pt,pp,series,predseries

def dist_calc(s1,s2):
    """
    Distance calculation for ABC runs.
    returns the mean of the RMS distances between series.
    """
    dists={}
    for n in s1.dtype.names:
        if n in s2:
            dists[n] = sqrt(mean((s1[n]-s2[n])**2))
    
    return mean(dists.values())

def main_loop(pars,ws=7,nc=20):
    """
    Run the full Analysis
        -pars: parameters for the
        -ws: Window size
        -nc: Number of inference cycles to run
    """
    K = 3000
    L = 1000
    global inits,tf
    hst = RTplot()
    hsp = RTplot()
    
    tf = ws*nc
    subnot = stats.uniform(.9,.1).rvs((tf,6))
    beta,alpha,sigma =pars
    duration = ws #size of the fitting window
    y = runModel(*[beta,alpha,sigma])
    #y *=subnot #apply noise due to sub-notification
    start = time.time()
    
    prior = {'theta':[],'phi':[]}
    for w in range(nc):
        print '==> Week # %s!'%w
        tf = duration
        if w > 0:
            inits = y[w*duration]
        pt,pp,series,predseries=doInference(K,L,duration,prior=prior,data=y[w*ws:w*ws+ws], method='SIR')
        f = open('week_%s'%w,'w')
        d = {'date':range(1,y.shape[0]+1),'S':y[:,0],'E':y[:,1],'I':y[:,2],'A':y[:,3],'R':y[:,4],'Ia':y[:,5]}
        cPickle.dump((pt,pp,series, d,predseries),f)
        f.close()
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


main_loop(pars,7,20)
plotWeekly(data=False,ws=7,nc=20)
#main_loop(pars,140,1)
#plotWeekly(data=False,ws=140,nc=1)
