# To change this template, choose Tools | Templates
# and open the template in the editor.
"""
Module with specialized plotting functions for the Melding results
"""

__author__="fccoelho"
__date__ ="$06/01/2010 11:24:20$"
__docformat__ = "restructuredtext en"

from itertools import cycle
from scipy import stats
import glob
import datetime
from numpy import *
import matplotlib.pyplot as P
from scipy.stats import gaussian_kde


def plot_series(tim,obs, series, names=[],title='series'):
    c = cycle(['b','g','r','c','m','y','k'])
    ser2={}
    for n in series[0].dtype.names:
        ser2[n] = concatenate([s[n] for s in series],axis=1)
    if not names:
        names = series[0].dtype.names
    for i,n in enumerate(names):
        #print n
        #P.subplot(l,1,i)
        co = c.next()
        P.plot(tim,[stats.scoreatpercentile(t,5) for t in ser2[n].T],co+'-.')
        P.plot(tim,[stats.scoreatpercentile(t,95) for t in ser2[n].T],co+'-.')
        P.plot(tim,[stats.scoreatpercentile(t,50) for t in ser2[n].T],co+'-',lw=2,label=n)
        #P.plot(tim,y[:,list(series.dtype.names).index(n)],'*',label=n+' obs.')
        if n in obs:
            P.plot(tim,obs[n],'*',label=n+' obs.')
    P.savefig(title+'.png')

def plot_pred(tim,series,y, fig,names=[],title='series'):
    c = cycle(['b','g','r','c','m','y','k'])
    if not names:
        names = series.dtype.names
    l = len(names)
    ax = fig.add_subplot(111)
    for i,n in enumerate(names):
        co = c.next()
        for b in [2.5,25,50,75,97.5]:
            if b == 50:
                st = 'k-'
            else:
                st = '%s:'%co
            try:
                lower = [stats.scoreatpercentile(t,0+b) for t in series[n].T]
                upper =[stats.scoreatpercentile(t,100-b) for t in series[n].T]
                ax.plot(tim,lower,st,tim,upper,st,label=n)
                ax.fill_between(tim,lower,upper,facecolor=co,alpha=0.12)#
            except:
                pass
        if n in y:
            #print y[n].shape, tim.shape
            try: #in the last window we run out of date points
                P.plot(tim,y[n],'^',label=n+' obs.')
            except:
                pass
    P.savefig(title+'.png')

def pred_new_cases(obs,series,weeks,names=[],ws=7):
    """
    Predicted total new cases in a window vs oserved.
    """
    fig =P.figure()
    P.title('Total new cases per window: predicted vs observed')
    if not names:
        names = series[0].dtype.names
    ax = fig.add_subplot(111)
    c = cycle(['b','g','r','c','m','y','k'])
    for n in names:
        if n in obs:
            co = c.next()
            ax.boxplot([sum(s[n],axis=1) for s in series] ,positions = range(weeks),notch=1,vert=1)
            ax.plot(range(weeks),[mean(sum(s[n],axis=1)) for s in series],'%s^'%co, label="Mean pred. %s"%n)
            ax.plot(range(weeks-1),[sum(obs[n][(w+1)*ws:(w+1)*ws+ws]) for w in range(weeks-1)],'%s-o'%co, label="obs. Prev")
    P.xlabel('windows')
    ax.legend(loc=0)

def plot_series2(tim,obs,series,names=[],title='Simulated vs Observed series',wl=7,lag=False):
    ser2={}
    for n in series[0].dtype.names:
        ser2[n] = concatenate([s[n] for s in series],axis=1)
    ls = ser2[n].shape[1]
    tim = tim[:ls]
    #print type (series)#.I.shape
    fig =P.figure()
    
    if not names:
        names = series[0].dtype.names
    c = cycle(['b','g','r','c','m','y','k'])
    if isinstance(tim[0], datetime.date):
        lag = datetime.timedelta(int(lag)*wl)
    else:
        lag = int(lag)*wl
    for i, n in enumerate(names):
        ax = fig.add_subplot(len(names), 1, i+1)
        co = c.next()
        if n in obs:
            ax.plot(tim,obs[n][:len(tim)],'%s+'%co, label="Observed %s"%n)
            #print len(tim),  ls
        ax.plot(array(tim)+lag,median(ser2[n],axis=0),'k-', label=n)
        lower = [stats.scoreatpercentile(t,2) for t in ser2[n].T]
        upper =[stats.scoreatpercentile(t,98) for t in ser2[n].T]
        ax.fill_between(array(tim)+lag,lower,upper,facecolor=co,alpha=0.6)
        if i < (len(names)-1):ax.xaxis.set_ticklabels([])
        ax.legend()
        if i == 0:
            ax.set_title(title)
    #ax.xaxis.set_visible(True)
    #P.title(title)
    P.xlabel('days')
    

def plot_par_series(tim,ptlist):
    P.figure()
    P.title('Parameters temporal variation')
    sq = sqrt(len(ptlist[0].dtype.names))
    r= floor(sq);c=ceil(sq)
    for i,n in enumerate(ptlist[0].dtype.names):
        P.subplot(r,c,i+1)
        P.boxplot([s[n] for s in ptlist],notch=1,positions=tim,vert=1)
        #P.errorbar(tim,[median(t[n]) for t in ptlist],yerr=[std(t[n]) for t in ptlist],label=n)
        P.ylabel(n)
    P.xlabel('Weeks')

def plot_par_violin(tim,ptlist, priors={}, bp=True):
    fig = P.figure()
    #P.title('Parameters temporal variation')
    sq = sqrt(len(ptlist[0].dtype.names))
    r= floor(sq);c=ceil(sq)
    if len(ptlist[0].dtype.names) == 3:
        r = 1; c = 3
    if priors:
        if len(tim)==1:
            tim  = [-1, 0]
        else:
            tim = [tim[0]-(tim[1]-tim[0])]+tim
    for i,n in enumerate(ptlist[0].dtype.names):
        ax = fig.add_subplot(r,c,i+1)
        violin_plot(ax,[priors[n]]+[s[n] for s in ptlist],tim,bp, True)
        P.ylabel(n)
    P.xlabel('Weeks')

def violin_plot(ax,data,positions,bp=False, prior = False):
    '''
    Create violin plots on an axis
    
    :Parameters:
        - `ax`: A subplot object
        - `data`: A list of data sets to plot
        - `positions`: x values to position the violins
        - `bp`: Whether to plot the boxplot on top.
        - `prior`: whether the first element of data is a Prior distribution.
    '''
    dist = max(positions)-min(positions)
    w = min(0.15*max(dist,1.0),0.5)
    i = 0
    for d,p  in zip(data, positions):
        if prior and i == 0:
            color = 'g'
        else:
            color = 'y'
        k = gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = arange(m,M,(M-m)/100.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        ax.fill_betweenx(x,p,v+p,facecolor=color,alpha=0.3)
        ax.fill_betweenx(x,p,-v+p,facecolor=color,alpha=0.3)
        i+=1
    if bp:
        ax.boxplot(data,notch=1,positions=positions,vert=1)

if __name__ == "__main__":
    print "Hello World";
