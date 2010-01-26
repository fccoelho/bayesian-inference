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
        co = c.next()
        ax.boxplot([sum(s[n],axis=1) for s in series] ,positions = range(weeks),notch=1,vert=1)
        ax.plot(range(weeks),[mean(sum(s[n],axis=1)) for s in series],'%s^'%co, label="Mean pred. %s"%n)
        ax.plot(range(weeks-1),[sum(obs[n][(w+1)*ws:(w+1)*ws+ws]) for w in range(weeks-1)],'%s-o'%co, label="obs. Prev")
    P.xlabel('weeks')
    ax.legend(loc=0)

def plot_series2(tim,obs,series,names=[],title='Simulated vs Observed series',ws=7,lag=False):
    ser2={}
    for n in series[0].dtype.names:
        ser2[n] = concatenate([s[n] for s in series],axis=1)
    #print type (series)#.I.shape
    fig =P.figure()
    P.title(title)
    if not names:
        names = series[0].dtype.names
    ax = fig.add_subplot(111)
    c = cycle(['b','g','r','c','m','y','k'])
    print len(tim), len(obs['I']),len(median(ser2['I'],axis=0))
    for n in names:
        co = c.next()
        if n in obs:
            ax.plot(tim,obs[n][:len(tim)],'%s+'%co, label="Observed %s"%n)
        ax.plot(array(tim)+ws*int(lag),median(ser2[n],axis=0),'k-')
        lower = [stats.scoreatpercentile(t,2) for t in ser2[n].T]
        upper =[stats.scoreatpercentile(t,98) for t in ser2[n].T]
        ax.fill_between(array(tim)+ws*int(lag),lower,upper,facecolor=co,alpha=0.6)
    P.xlabel('days')
    ax.legend()

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

def plot_par_violin(tim,ptlist):
    fig = P.figure()
    #P.title('Parameters temporal variation')
    sq = sqrt(len(ptlist[0].dtype.names))
    r= floor(sq);c=ceil(sq)
    for i,n in enumerate(ptlist[0].dtype.names):
        ax = fig.add_subplot(r,c,i+1)
        violin_plot(ax,[s[n] for s in ptlist],tim,bp=True)
        P.ylabel(n)
    P.xlabel('Weeks')

def violin_plot(ax,data,positions,bp=False):
    '''
    Create violin plots on an axis
    '''
    dist = max(positions)-min(positions)
    w = min(0.15*max(dist,1.0),0.5)
    for d,p in zip(data,positions):
        k = gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = arange(m,M,(M-m)/100.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
        ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
    if bp:
        ax.boxplot(data,notch=1,positions=positions,vert=1)

if __name__ == "__main__":
    print "Hello World";
