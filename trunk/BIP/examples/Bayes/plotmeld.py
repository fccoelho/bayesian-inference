"""
Plot the results of the recursive parameter estimation simulations
and resulting time-series
"""

#
# Copyright 2009- by Flávio Codeço Coelho
# License gpl v3
#
import cPickle
from itertools import cycle
from scipy import stats
import glob
from numpy import *
import matplotlib.pyplot as P
import Gnuplot
from scipy.stats import gaussian_kde




        
def plotRaHist(arr, title=''):
    '''
    Plots a record array
    as a panel of histograms
    '''
    nv = len(arr.dtype.names)
    fs = (ceil(sqrt(nv)),floor(sqrt(nv))+1) #figure size
    P.figure()
    for i,n in enumerate(arr.dtype.names):
        P.subplot(nv/2+1,2,i+1)
        P.hist(arr[n],bins=50, normed=True, label=n)
        P.legend()
    P.savefig(title+'.png')


def plotParSeries(tim,ptlist):
    P.figure()
    P.title('Parameters temporal variation')
    for i,n in enumerate(ptlist[0].dtype.names):
        P.subplot(3,1,i+1)
        P.boxplot([s[n] for s in ptlist],notch=1,positions=tim,vert=1)
        #P.errorbar(tim,[median(t[n]) for t in ptlist],yerr=[std(t[n]) for t in ptlist],label=n)
        P.ylabel(n)
    P.xlabel('Weeks')
    
def plot_par_violin(tim,ptlist):
    fig = P.figure()
    #P.title('Parameters temporal variation')
    for i,n in enumerate(ptlist[0].dtype.names):
        ax = fig.add_subplot(3,1,i+1)
        violin_plot(ax,[s[n] for s in ptlist],tim,bp=True)
        P.ylabel(n)
    P.xlabel('Weeks')
    
def violin_plot(ax,data,positions,bp=False):
    '''
    create violin plots on an axis
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

def plotSeries(tim,series,y, names=[],title='series'):
    c = cycle(['b','g','r','c','m','y','k'])
    if not names:
        names = series.dtype.names
    l = len(names)
    for i,n in enumerate(names):
        #P.subplot(l,1,i)
        co = c.next()
        P.plot(tim,[stats.scoreatpercentile(t,5) for t in series[n].T],co+'-.')
        P.plot(tim,[stats.scoreatpercentile(t,95) for t in series[n].T],co+'-.')
        P.plot(tim,[stats.scoreatpercentile(t,50) for t in series[n].T],co+'-',lw=2,label=n)
        #P.plot(tim,y[:,list(series.dtype.names).index(n)],'*',label=n+' obs.')
        if n in y:
            P.plot(tim,y[n],'*',label=n+' obs.')
    P.savefig(title+'.png')

def plotPred(tim,series,y, fig,names=[],title='series'):
    c = cycle(['b','g','r','c','m','y','k'])
    if not names:
        names = series.dtype.names
    l = len(names)
    ax = fig.add_subplot(111)
    for i,n in enumerate(names):
        co = c.next()
        for b in range(10,90,10):
            if b == 50:
                st = 'k-'
            else:
                st = '%s:'%co
            lower = [stats.scoreatpercentile(t,0+b) for t in series[n].T]
            upper =[stats.scoreatpercentile(t,100-b) for t in series[n].T]
            ax.plot(tim,lower,st,tim,upper,st,label=n)
            ax.fill_between(tim,lower,upper,facecolor=co,alpha=0.12)#
        if n in y:
            P.plot(tim,y[n],'^',label=n+' obs.')
    P.savefig(title+'.png')

def predNewCases(obs,series,weeks,fig, names=[],ws=7):
    """
    Predicted total new cases in a week vs oserved.
    """
    fig =P.figure()
    P.title('Total new cases per week: predicted vs observed')
    if not names:
        names = series.dtype.names
    ax = fig.add_subplot(111)
    
    for n in names:
        ax.boxplot([sum(s[n],axis=1) for s in series] ,positions = range(weeks),notch=1,vert=1)
        ax.plot(range(weeks),[mean(sum(s[n],axis=1)) for s in series],'^', label="Mean pred. %s"%n)
        ax.plot(range(weeks-1),[sum(obs[n][(w+1)*ws:(w+1)*ws+ws]) for w in range(weeks-1)],'r-o', label="obs. Prev")
    P.xlabel('weeks')
    ax.legend(loc=0)

def plot_series2(tim,obs,series,names=[],tit='Simulated vs Observed series',ws=7,lag=False):
    ser2={}
    for n in series[0].dtype.names:
        ser2[n] = concatenate([s[n] for s in series],axis=1)
    #print type (series)#.I.shape
    fig =P.figure()
    P.title(tit)
    if not names:
        names = ser2.keys()
    ax = fig.add_subplot(111)
    c = cycle(['b','g','r','c','m','y','k'])
    print len(tim), len(obs['I']),len(median(ser2['I'],axis=0))
    for n in names:
        co = c.next()
        ax.plot(tim,obs[n],'%s+'%co, label="Observed %s"%n)
        ax.plot(array(tim)+ws*int(lag),median(ser2[n],axis=0),'k-')
        lower = [stats.scoreatpercentile(t,2) for t in ser2[n].T]
        upper =[stats.scoreatpercentile(t,98) for t in ser2[n].T]
        ax.fill_between(array(tim)+ws*int(lag),lower,upper,facecolor=co,alpha=0.6)
    P.xlabel('days')
    ax.legend()    
    
def plotWeekly(data=True,ws=7,nc=1):
    pt,pp,series,predseries,fi=[],[],[],[],[]
    if not data:
        plot_wo_data(pt,pp,series,predseries,fi,ws,nc)
        return
    for w in range(nc):
        fn = 'weekd_%s'%w
        print fn
        f = open(fn,'r')
        a,b,c,obs,pred, samples = cPickle.load(f)
        f.close()
    
        pt.append(a)
        pp.append(b)
        series.append(c)
        predseries.append(pred)
        fi.append(a.beta*b.I*b.S)
        #predictions(obs,pred,w+1)
    incid = concatenate([s.I for s in series], axis=1)
    #plotIncidence(incid)
    #plotParHistSeries(pt)
    #P.errorbar(range(len(pt)),[mean(i) for i in fi],yerr=[std(i) for i in fi],label="force of infection")
    plotParSeries(range(len(pt)),pt)
    plot_par_violin(range(len(pt)),pt)
    P.figure()
    for w in range(nc):
        ow = dict([(k,v[w*ws:w*ws+ws]) for k,v in obs.items()])
        plotSeries(ow['date'],series[w], ow,['I','Ia','A','R'])
        if w==0:
            P.legend()
            P.title('Weekly Adjusted Series' )
    fig=P.figure()
    for w in range(nc-1):
        owp = dict([(k,v[(w+1)*ws:(w+1)*ws+ws]) for k,v in obs.items()])
        #print y.shape, predseries.shape
        plotPred(owp['date'],predseries[w], owp,fig,['I','Ia'])
        
        if w==0:
            #P.legend()
            P.title('Weekly Predicted Series' )
    predNewCases(obs,predseries,nc,fig,['I','Ia'],ws)
    P.show()

def plot_wo_data(pt,pp,series,predseries,fi,ws=7,nc=1):
    for w in range(nc):
        fn = 'week_%s'%w
        print "Reading ",fn
        f = open(fn,'r')
        a,b,c,obs,pred = cPickle.load(f)
        f.close()
        
        pt.append(a)
        pp.append(b)
        series.append(c)
        predseries.append(pred)
        fi.append(a[r'$\beta$']*b.I*b.S)
    #predictions(obs,pred,w+1)
    incid = concatenate([s.I for s in series], axis=1)
    #plotIncidence(incid)
    #plotParHistSeries(pt)
    #P.errorbar(range(len(pt)),[mean(i) for i in fi],yerr=[std(i) for i in fi],label="force of infection")
    plotParSeries(range(len(pt)),pt)
    plot_par_violin(range(len(pt)),pt)
    P.figure()
    for w in range(nc):
        ow = dict([(k,v[w*ws:w*ws+ws]) for k,v in obs.items()])
        plotSeries(ow['date'],series[w], ow,['I','A'])
        if w==0:
            P.legend()
            P.title('Weekly Adjusted Series' )
    fig=P.figure()
    for w in range(nc-1):
        owp = dict([(k,v[(w+1)*ws:(w+1)*ws+ws]) for k,v in obs.items()])
        #print y
        plotPred(owp['date'],predseries[w], owp,fig,['I','A'])
        
        if w==0:
            #P.legend()
            P.title('Weekly Predicted Series' )
    predNewCases(obs,predseries,nc,fig,['I','A'],ws)
    plot_series2(range(nc*ws),obs,series,names=['I','A'])
    plot_series2(range(nc*ws),obs,predseries,names=['I','A'],tit='Predicted vs Observed series',lag=True)
    P.show()

if __name__=="__main__":
    plotWeekly(0,7,20)
    #plotWeekly()
