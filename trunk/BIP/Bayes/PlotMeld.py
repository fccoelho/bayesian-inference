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
from matplotlib.dates import date2num
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import gaussian_kde
import pdb



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
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
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

def pred_new_cases(obs,series,weeks,names=[], title='Total new cases per window: predicted vs observed' ,ws=7):
    """
    Predicted total new cases in a window vs oserved.
    """
    fig =P.gcf()
    P.title(title)
    if not names:
        names = series[0].dtype.names
    ax = P.gca()#fig.add_subplot(111)
    c = cycle(['b','g','r','c','m','y','k'])
    if 'time' in obs: #setting the xlabel 
        x = date2num([obs['time'][ws*i] for i in range(1, weeks)])
        ax.xaxis_date()
    else:
        x = arange(1, weeks)
    sc= 1 if len(series) ==1 else 5
    W = min(0.5*max(len(x),1.0),0.5)*sc
    for n in names:
        if n in obs:
            co = c.next()
            print len(x),  len([mean(sum(s[n],axis=1)) for s in series]),  type(x)
            ax.plot([x[7]]+x.tolist(), [mean(sum(s[n],axis=1)) for s in series],'%s^'%co, label="Mean pred. %s"%n)
            ax.plot(x,[nansum(obs[n][(w+1)*ws:(w+1)*ws+ws]) for w in range(weeks-1)],'%s-o'%co, label="obs. Prev")
            ax.boxplot([nansum(s[n],axis=1) for s in series] ,positions = x, widths=W,notch=1,vert=1)
    #P.xlabel('windows')
    #ax.legend(loc=0)
    if 'time' in obs: 
        fig.autofmt_xdate()

def plot_series2(tim,obs,series,names=[],title='Simulated vs Observed series',wl=7,lag=False):
    ser2={}
    for n in series[0].dtype.names:
        ser2[n] = concatenate([s[n] for s in series],axis=1)
    ls = ser2[n].shape[1]
    tim = tim[:ls]
    #print type (series)#.I.shape
    fig =P.gcf()
    if not names:
        names = series[0].dtype.names
    c = cycle(['b','g','r','c','m','y','k'])
    if isinstance(tim[0], datetime.date):
        lag = datetime.timedelta(int(lag)*wl)
    else:
        lag = int(lag)*wl
    for i, n in enumerate(names):
        ax = fig.add_subplot(len(names), 1, i+1)
        ax.grid(True)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        P.setp(ax.get_yticklabels(),fontsize=8)
        P.setp(ax.get_xticklabels(),fontsize=8)
        if isinstance(tim[0], datetime.date):
            ax.xaxis_date()
        co = c.next()
        if n in obs:
            ax.plot(tim,obs[n][:len(tim)],'o', label=r"$Observed\; %s$"%n)
        #pdb.set_trace()
        ax.plot(array(tim)+lag,median(ser2[n],axis=0),'k-', label=r"$median\; %s$"%n)
        ax.plot(array(tim)+lag,mean(ser2[n],axis=0),'k--', label=r"$mean\; %s$"%n)
        lower = [stats.scoreatpercentile(t,2.5) for t in ser2[n].T]
        upper =[stats.scoreatpercentile(t,97.5) for t in ser2[n].T]
        if len(series)>1: #in the case of iterative simulations
            dif = (array(upper)-array(lower))
            dif  = dif/max(dif)*10
            pe, va = peakdet(dif, 1)
            xp = [0]+ pe[:, 0].tolist()+[len(lower)-1]
            lower = interp(range(len(lower)), xp, array(lower)[xp]) # valley-to-valley interpolated band
            upper = interp(range(len(upper)), xp, array(upper)[xp])#peak-to-peak interpolated band
        ax.fill_between(array(tim)+lag,lower,upper,facecolor=co,alpha=0.2)
        #ax.fill_between(array(tim)+lag,lower,upper,facecolor='k',alpha=0.1)
        if i < (len(names)-1):ax.xaxis.set_ticklabels([])
        ax.legend()
        if i == 0:
            ax.set_title(title)
    #ax.xaxis.set_visible(True)
    #P.title(title)
    P.xlabel('days')
    if isinstance(tim[0], datetime.date):
        fig.autofmt_xdate()
    

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
    P.xlabel('Windows')

def plot_par_violin(tim,ptlist, priors={}, bp=True):
    fig = P.figure()
    #P.title('Parameters temporal variation')
    sq = sqrt(len(ptlist[0].dtype.names))
    ad = 1 if sq%1 >0.5 else 0
    r= floor(sq)+ad;c=ceil(sq)
    if len(ptlist[0].dtype.names) == 3:
        r = 3; c = 1
    if priors:
        if isinstance(tim[0], datetime.date):
            pdate = tim[0]-datetime.timedelta(1) if len(tim)==1 else tim[0]-(tim[1]-tim[0])
            tim =[pdate]+tim.tolist()
        else:
            if len(tim)==1:
                tim  = [-1, 0]
            else:
                tim = [tim[0]-(tim[1]-tim[0])]+tim
    for i,n in enumerate(ptlist[0].dtype.names):
        ax = fig.add_subplot(r,c,i+1)
        ax.grid(True)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        P.setp(ax.get_yticklabels(),fontsize=8)
        P.setp(ax.get_xticklabels(),fontsize=8)
        violin_plot(ax,[priors[n]]+[s[n] for s in ptlist],tim,bp, True)
        P.ylabel(n)
    #P.xlabel('Windows')
    if isinstance(tim[0], datetime.date):
        fig.autofmt_xdate()

def violin_plot(ax,data,positions,bp=False, prior = False):
    '''
    Create violin plots on an axis
    
    :Parameters:
        - `ax`: A subplot object
        - `data`: A list of data sets to plot
        - `positions`: x values to position the violins. Can be datetime.date objects.
        - `bp`: Whether to plot the boxplot on top.
        - `prior`: whether the first element of data is a Prior distribution.
    '''
    sc = 1
    dist = len(positions)
    if isinstance(positions[0], datetime.date):
        ax.xaxis_date()
        positions = date2num(positions)
        sc = 5 if (dist>2 ) else 1
        #print sc
    w = min(0.5*max(dist,1.0),0.5)*sc
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
        ax.boxplot(data,notch=1,positions=positions,widths=w, vert=1)

def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    Currently returns two lists of tuples, but maybe arrays would be better
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    % [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    % maxima and minima ("peaks") in the vector V.
    % MAXTAB and MINTAB consists of two columns. Column 1
    % contains indices in V, and column 2 the found values.
    %
    % With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    % in MAXTAB and MINTAB are replaced with the corresponding
    % X-values.
    %
    % A point is considered a maximum peak if it has the maximal
    % value, and was preceded (to the left) by a value lower by
    % DELTA.
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
 
    return array(maxtab), array(mintab)

if __name__ == "__main__":
    print "Hello World";
