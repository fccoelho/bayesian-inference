# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
# Name:        Melding.py
# Purpose:     The Bayesian melding Class provides
#                   uncertainty analyses for deterministic models.
#
# Author:      Flávio Codeço Coelho
#
# Created:     2003/08/10
# Copyright:   (c) 2003-2008 by the Author
# Licence:     GPL
#-----------------------------------------------------------------------------
from numpy.core.records import recarray
import psyco
psyco.full()
import sys, os
import cPickle as CP
import like
import pylab as P
import scipy.stats.kde as kde
from scipy import stats
import numpy
from numpy import *
from time import time
from numpy.random import normal, randint,  random,  uniform 
import lhs
if sys.version.startswith('2.5'):
    from processing import Pool
else:
    from multiprocessing import Pool

__docformat__ = "restructuredtext en"


class Meld:
    """
    Bayesian Melding class
    """
    def __init__(self,  K,  L, model, ntheta, nphi, alpha = 0.5, verbose = False ):
        """
        Initializes the Melding class.
        
        :Parameters:
            - `K`: Number of replicates of the model run. Also determines the prior sample size.
            - `L`: Number of samples from the Posterior distributions. Usually 10% of K.
            - `model`: Callable taking theta as argument and returning phi = M(theta).
            - `ntheta`: Number of inputs to the model (parameters).
            - `nphi`: Number of outputs of the model (State-variables)
        """
        self.K = K
        self.L = L
        self.verbose = verbose
        self.model = model
        self.likelist = [] #list of likelihoods
        self.q1theta = recarray(K,formats=['f8']*ntheta) #Theta Priors (record array)
        self.post_theta = recarray(L,formats=['f8']*ntheta) #Theta Posteriors (record array)
        self.q2phi = recarray(K,formats=['f8']*nphi) #Phi Priors (record array)
        self.phi = recarray(K,formats=['f8']*nphi) #Phi model-induced Priors (record array)
        self.q2type = [] #list of distribution types
        self.post_phi = recarray(L,formats=['f8']*nphi) #Phi Posteriors (record array)
        self.ntheta = ntheta
        self.nphi = nphi
        self.alpha = alpha #pooling weight of user-provided phi priors
        self.done_running = False
#        self.po = Pool() #pool of processes for parallel processing
    
    def setPhi(self, names, dists=[stats.norm], pars=[(0, 1)], limits=[(-5,5)]):
        """
        Setup the models Outputs, or Phi, and generate the samples from prior distributions 
        needed for the melding replicates.
        
        :Parameters:
            - `names`: list of string with the names of the variables.
            - `dists`: is a list of RNG from scipy.stats
            - `pars`: is a list of tuples of variables for each prior distribution, respectively.
            - `limits`: lower and upper limits on the support of variables.
        """
        if len(names) != self.nphi:
            raise ValueError("Number of names(%s) does not match the number of output variables(%s)."%(len(names),self.nphi))
        self.q2phi.dtype.names = names
        self.phi.dtype.names = names
        self.post_phi.dtype.names = names
        self.limits = limits
        for n,d,p in zip(names,dists,pars):
            self.q2phi[n] = lhs.lhs(d,p,self.K)
            self.q2type.append(d.name)


        
    def setTheta(self, names, dists=[stats.norm], pars=[(0, 1)]):
        """
        Setup the models inputs and generate the samples from prior distributions 
        needed for the dists the melding replicates.
        
        :Parameters:
            - `names`: list of string with the names of the parameters.
            - `dists`: is a list of RNG from scipy.stats
            - `pars`: is a list of tuples of parameters for each prior distribution, respectivelydists
        """
        self.q1theta.dtype.names = names
        self.post_theta.dtype.names = names
        if os.path.exists('q1theta'):
            self.q1theta = CP.load(open('q1theta','r'))
        else:
            for n,d,p in zip(names,dists,pars):
                self.q1theta[n] = lhs.lhs(d,p,self.K)
        
    def setThetaFromData(self,names,data):
        """
        Setup the model inputs and set the prior distributions from the vectors
        in data.
        This method is to be used when the prior distributions are available in 
        the form of a sample from an empirical distribution such as a bayesian
        posterior.
        In order to expand the samples provided, K samples are generated from a
        kernel density estimate of the original sample.
        
        :Parameters:
            - `names`: list of string with the names of the parameters.
            - `data`: list of vectors. Samples of a proposed distribution
        """
        self.q1theta.dtype.names = names
        self.post_theta.dtype.names = names
        if os.path.exists('q1theta'):
            self.q1theta = CP.load(open('q1theta','r'))
        else:
            for n,d in zip(names,data):
                self.q1theta[n] = kde.gaussian_kde(d).resample(self.K)

    def setPhiFromData(self,names,data,limits):
        """
        Setup the model outputs and set their prior distributions from the
        vectors in data.
        This method is to be used when the prior distributions are available in
        the form of a sample from an empirical distribution such as a bayesian
        posterior.
        In order to expand the samples provided, K samples are generated from a
        kernel density estimate of the original sample.

        :Parameters:
            - `names`: list of string with the names of the variables.
            - `data`: list of vectors. Samples of the proposed distribution.
            - `limits`: list of tuples (ll,ul),lower and upper limits on the support of variables.
        """
        self.q2phi.dtype.names = names
        self.phi.dtype.names = names
        self.post_phi.dtype.names = names
        self.limits = limits
        for n,d in zip(names,data):
            self.q2phi[n] = kde.gaussian_kde(d).resample(self.K)
            self.q2type.append('empirical')

    def addData(self, data, model, limits,l=1024, **kwargs):
        """
        Calculates the likelihood functions of the dataset presented and add to 
        self.likelist
        Likelihood function is a vector of lenght l
        
        :Parameters:
            -  `data`: vector containing observations on a given variable.
            -  `model`: string with the name of the distribution of the variable
            -  `limits`: (ll,ul) tuple with lower and upper limits for the variable
            -  `l`: Length (resolution) of the likelihood vector
        """
        n = len(data) # Number of data points
        data = array(data)
        (ll,ul) = limits #limits for the parameter space
        step = (ul-ll)/float(l)
        
        if model == 'normal': # In this case, L is a function of the mean. SD is set to the SD(data)
            sd = std(data) #standard deviation of data
            prec = 1/sd #precision of the data
            res = array([exp(like.Normal(data,mu,prec)) for mu in arange(ll,ul,step)])  
            lik = res/max(res) # Likelihood function   
            print max(lik), min(lik)
        elif model == 'exponential':
            res = [lamb**n*exp(-lamb*sum(data)) for lamb in arange(ll,ul,step)]
            lik = array(res)/max(array(res))
        elif model == 'beta':
            # TODO: Make sure pars is passed as an extra parameter
            res = [exp(like.Beta(data,*kwargs['pars'])) for i in arange(ll,ul,step)]
            lik = array(res)/max(array(res))
        elif model == 'bernoulli':
            if ll<0 or ul>1:
                print "Parameter p of the bernoulli is out of range[0,1]"
            res = [exp(like.Bernoulli(data,p)) for p in arange(ll,ul,step)]
            lik = array(res)/max(array(res))
            
        elif model == 'poisson':
            res = [exp(like.Poisson(data,lb)) for lb in arange(ll,ul,step)]
            lik = array(res)/max(array(res))
        
        elif model == 'lognormal':
            sd = std(data) #standard deviation of data
            prec = 1/sd #precision of the data
            res = [exp(like.Lognormal(data,mu,prec)) for mu in arange(ll,ul,step)]
            lik = array(res)/max(array(res))    
        else:
            print 'Invalid distribution type. Valid distributions: normal,lognormal, exponential, bernoulli and poisson'
        self.likelist.append(lik)
        return lik
        
    def run(self,*args):
        """
        Runs the model through the Melding inference.model
        model is a callable which return the output of the deterministic model,
        i.e. the model itself.
        The model is run self.K times to obtain phi = M(theta).
        """
        
        for i in xrange(self.K):
            theta = [self.q1theta[n][i] for n in self.q1theta.dtype.names]
            r = self.po.applyAsync(self.model, theta)
            self.phi[i]= r.get()[-1]#self.model(*theta)[-1] #phi is the last point in the simulation

        self.done_running = True
        
    def getPosteriors(self,t=1):
        """
        Updates the the posteriors of the model for the last time step.
        Returns two record arrays:
        - The posteriors of the Theta
        - the posterior of Phi last t values of time-series. self.L by `t` arrays.

        :Parameters:
            - `t`: length of the time-series to return as posterior.
        """
        if not self.done_running:
            return
        if t > 1:
            self.post_phi = recarray((self.L,t),formats=['f8']*self.nphi)
            self.post_phi.dtype.names = self.phi.dtype.names
        #random indices for the marginal posteriors of theta
        po = Pool()
        #pti = randint(0,self.L,size=(self.ntheta,self.L))
        pti = lhs.lhs(stats.randint,(0,self.L),siz=(self.ntheta,self.L))
        for i in xrange(self.L):#Monte Carlo with values of the posterior of Theta
            r = po.applyAsync(self.model,[self.post_theta[n][pti[j,i]] for j,n in enumerate(self.post_theta.dtype.names)] )
            if t == 1:
                self.post_phi[i] = r.get()[-1]
            else:
                self.post_phi[i]= [tuple(l) for l in r.get()[-t:]]
            if not i%100:
                print "==> L = %s"%i

        po.close()
        po.join()
        return self.post_theta, self.post_phi

    def filtM(self,cond,x,limits):
        '''
        Multiple condition filtering.
        Remove values in x[i], if corresponding values in
        cond[i] are less than limits[i][0] or greater than
        limits[i][1].

        :Parameters:
            - `cond`: is an array of conditions.
            - `limits`: is a list of tuples (ll,ul) with length equal to number of lines in `cond` and `x`.
            - `x`: array to be filtered.
        '''
        cond = array(cond)
        cnd = ones(cond.shape[1],int)
        for i,j in zip(cond,limits):
            ll = j[0]
            ul = j[1]
            #print cond.shape,cnd.shape,i.shape,ll,ul
            cnd = cnd & less(i,ul) & greater(i,ll)
        f = compress(cnd,x, axis=1)
        return f

    def basicfit(self,s1,s2):
        '''
        Calculates a basic fitness calculation between a model-
        generated time series and a observed time series.
        it uses a normalized RMS variation.

        :Parameters:
            - `s1`: model-generated time series. record array.
            - `s2`: observed time series. dictionary with keys matching names of s1
        '''
        fit = []
        for k in s2.keys():
            if s2[k] == [] or (not s2[k].any()):
                continue #no observations for this variable
            e = sqrt(mean((s1[k]-s2[k])**2.))
            fit.append(e) #min to guarantee error is bounded to (0,1)

        return mean(fit) #mean r-squared
        

    def logPooling(self,phi):
        """
        Returns the probability associated with each phi[i]
        on the pooled pdf of phi and q2phi.

        :Parameters:
            - `phi`: prior of Phi induced by the model and q1theta.
        """
        
#       Estimating the multivariate joint probability densities
        #print phi[phi.dtype.names[0]].shape
        phidens = stats.gaussian_kde(array([phi[n][:,-1] for n in phi.dtype.names]))
        q2dens = stats.gaussian_kde(array([self.q2phi[n] for n in self.q2phi.dtype.names]))
#       Determining the pooled probabilities for each phi[i]
#        qtilphi = zeros(self.K)
        lastp = array([list(phi[i,-1]) for i in xrange(self.K)])
#        print lastp,lastp.shape
        qtilphi = (phidens.evaluate(lastp.T)**(1-self.alpha))*q2dens.evaluate(lastp.T)**self.alpha
        return qtilphi/sum(qtilphi)

    def abcRun(self,fitfun=None, data={}, t=1,savetemp=False):
        """
        Runs the model for inference through Approximate Bayes Computation
        techniques. This method should be used as an alternative to the sir.
        
        :Parameters:
             - `fitfun`: Callable which will return the goodness of fit of the model to data as a number between 0-1, with 1 meaning perfect fit
             - `t`: number of time steps to retain at the end of the of the model run for fitting purposes.
             - `data`: dict containing observed time series (lists of length t) of the state variables. This dict must have as many items the number of state variables, with labels matching variables names. Unorbserved variables must have an empty list as value.
             - `savetemp`: Should temp results be saved. Useful for long runs. Alows for resuming the simulation from last sa
        """
        if not fitfun:
            fitfun = self.basicfit
        if savetemp:
            CP.dump(self.q1theta,open('q1theta','w'))
#       Running the model ==========================
        if os.path.exists('phi.temp'):
            phi,j = CP.load(open('phi.temp','r'))
        else:
            j=0
            phi = recarray((self.K,t),formats=['f8']*self.nphi, names = self.phi.dtype.names)
        for i in xrange(j,self.K):
            theta = [self.q1theta[n][i] for n in self.q1theta.dtype.names]
            r = self.po.applyAsync(self.model, theta)
            phi[i]= [tuple(l) for l in r.get()[-t:]]# #phi is the last t points in the simulation
            if i%100 == 0:
                print "==> K = %s"%i
                if savetemp:
                    CP.dump((phi,i),open('phi.temp','w'))
        if savetemp: #If all replicates are done, clear temporary save files.
            os.unlink('phi.temp')
            os.unlink('q1theta')

        print "==> Done Running the K replicates\n"
        qtilphi = self.logPooling(phi) #vector with probability of each phi[i] belonging to qtilphi
        qtilphi = nan_to_num(qtilphi)
        print 'max(qtilphi): ', max(qtilphi)
#      
#        calculate weights
        w = [fitfun(phi[i],data) for i in xrange(phi.shape[0])]
        w /=sum(w)
        w = 1-w
        print "w=",w, mean(w), var(w)
        print
        print 'qtilphi=',qtilphi
        # Resampling Thetas
        w = nan_to_num(w)
        w = array(w)*qtilphi
        w /=sum(w)
        w = nan_to_num(w)
        print 'max(w): ',max(w)
#        for n in phi.dtype.names:
#            P.plot(mean(phi[n],axis=0),label=n)
#        P.figure()
#        P.plot(w,label='w')
#        P.plot(qtilphi,label='qtilphi')
#        P.title('Resampling vector(w) and pooled prior on Phi')
#        P.legend()
        if sum(w) == 0.0:
            sys.exit('Resampling weights are all zero, please check your model or data.')
        j = 0
        while j < self.L: # Extract L samples from q1theta
            i=randint(0,w.size)# Random position of w and q1theta
            if random()<= w[i]:
                self.post_theta[j] = self.q1theta[i]# retain the sample according with resampling prob.
                j+=1
        

        self.done_running = True

    def sir(self, data={}, t=1,savetemp=False):
        """
        Run the model output through the Sampling-Importance-Resampling algorithm

        :Parameters:
            - `data`: observed time series on the model's output
            - `t`: length of the observed time series
            - `savetemp`: Boolean. create a temp file?
        """
        qtilphi,phi = self.runModel(savetemp,t)

#        Calculating the likelihood of each phi[i] considering the observed data
        tau = 0.1
        lik = zeros(self.K)
        t0=time()
        for i in xrange(self.K):
            l=1
            for n in data.keys():
                if isinstance(data[n],list) and data[n] == []: 
                    continue #no observations for this variable
                elif isinstance(data[n],numpy.ndarray) and (not data[n].any()):
                    continue #no observations for this variable
                p = phi[n]
#                po = Pool()
#                ps = [po.applyAsync(like.Normal,(data[n][m], j, tau)) for m,j in enumerate(p[i])]
#                l = product([p.get() for p in ps])
#                po.close()
#                po.join()
                l *= product([like.Normal(data[n][m], j, tau) for m,j in enumerate(p[i])])
            lik[i]=l
        print "==> Done Calculating Likelihoods (took %s seconds)"%(time()-t0)
#        Calculating the weights
        w = nan_to_num(qtilphi*lik)
        w = nan_to_num(w/sum(w))

        if sum(w) == 0.0:
            sys.exit('Resampling weights are all zero, please check your model or data.')
        j = 0
        while j < self.L: # Extract L samples from q1theta
            i=randint(0,w.size)# Random position of w and q1theta
            if random()<= w[i]:
                self.post_theta[j] = self.q1theta[i]# retain the sample according with resampling prob.
                j+=1
        self.done_running = True


    def runModel(self,savetemp,t=1):
        '''
        Handles running the model self.K times keeping a temporary savefile for 
        resuming calculation in case of interruption.

        :Parameters:
            - `savetemp`: Boolean. create a temp file?
        '''
        if savetemp:
            CP.dump(self.q1theta,open('q1theta','w'))
#       Running the model ==========================
        
                
        if os.path.exists('phi.temp'):
            phi,j = CP.load(open('phi.temp','r'))
        else:
            j=0
            phi = recarray((self.K,t),formats=['f8']*self.nphi, names = self.phi.dtype.names)
        def cb(r):
            '''
            callback function for the asynchronous model runs
            '''
            if t == 1:
                phi[r[1]] = (r[0][-1],)
            else:
                phi[r[1]] = [tuple(l) for l in r[0][-t:]]# #phi is the last t points in the simulation

        po = Pool()
        t0=time()
        for i in xrange(j,self.K):
            theta = [self.q1theta[n][i] for n in self.q1theta.dtype.names]
            r = po.applyAsync(enumRun,(self.model,theta,i),callback=cb)
#            r = po.applyAsync(self.model,theta)
#            if t == 1:
#                phi[i] = (r.get()[-1],)
#            else:
#                phi[i] = [tuple(l) for l in r.get()[-t:]]# #phi is the last t points in the simulation
            if i%100 == 0 and self.verbose:
                print "==> K = %s"%i
                if savetemp:
                    CP.dump((phi,i),open('phi.temp','w'))
        if savetemp: #If all replicates are done, clear temporary save files.
            os.unlink('phi.temp')
            os.unlink('q1theta')
        po.close()
        po.join()
        print "==> Done Running the K replicates (took %s seconds)\n"%(time()-t0)
        t0 = time()
        qtilphi = self.logPooling(phi) #vector with probability of each phi[i] belonging to qtilphi
        print "==> Done Running the Log Pooling (took %s seconds)\n"%(time()-t0)
        print qtilphi,'max(qtilphi): ', max(qtilphi)
        qtilphi = nan_to_num(qtilphi)
        return qtilphi,phi
def enumRun(model,theta,k):
    res =model(*theta)
    return (res,k)

def model(r, p0, n=1):
    """
    Model (r,p0, n=1)
    Simulates the Population dynamic Model (PDM) Pt = rP0
    for n time steps.
    P0 is the initial population size. 
    Example model for testing purposes.
    """
#    print "oi"
    Pt = zeros(n, float) # initialize the output vector
    P = p0
    for i in xrange(n):
        Pt[i] = r*P
        P = Pt[i]
    
    return Pt

def Run(k):
    """
    Run (k)
    Draw k samples of Theta from its prior distribution, run the model with it 
    and obtain phi = M(theta). For testing purposes only.
    """
    po = Pool()
#---q1theta---------------------------------------------------------------------
#---Priors for the theta (model parameters)--------------------
    r = lhs.lhs(stats.uniform, [2, 4], k)
    p0 = lhs.lhs(stats.uniform,[0,5],k)
    q1theta = (r, p0)
#-------------------------------------------------------------------------------
    phi=zeros(k, float)
    #print r.shape, p0.shape
    for i in xrange(k):
        re = po.applyAsync(model,(r[i], p0[i]))
        phi[i] = re.get()[-1]#model(r[i], p0[i])[-1] # Sets phi[i] to the last point of the simulation
        
    
    return phi, q1theta

def KDE(x, (ll, ul)=('',''),res=1024.):
    """
    KDE(x)
    performs a kernel density estimate using the scipy gaussian density
    if (ll,ul), enforce limits for the distribution's support.
    Returns a dictionary.
    """
    #r.assign("x", x)
    
    if ll :
        rn=arange(ll,ul,(ul-ll)/res)
        #print x.shape,rn.shape
        est = kde.gaussian_kde(x.ravel()).evaluate(rn)
        #r.assign("ll", ll)
        #r.assign("ul", ul)
        #est = r('density(x,from=ll, to=ul)') #trims the density borders
    else:
        ll = min(x)
        ul = max(x)
        rn=arange(ll,ul,(ul-ll)/res)
        est = kde.gaussian_kde(x).evaluate(rn)
        #est = r('density(x)')
        print 'No - KDE'
    return {'y':est,'x':rn}


def Likeli(data, dist, limits,**kwargs):
    """
    Generates the likelihood function of data given dist.
    limits is a tuple setting the interval of the parameter space that will
    be used as the support for the Likelihood function.
    returns a vector (1024 elements).
    """
    n = len(data) # Number of data points
    data = array(data)
    (ll,ul) = limits #limits for the parameter space
    step = (ul-ll)/1024.
    
    if dist == 'normal': # In this case, L is a function of the mean. SD is set to the SD(data)
        sd = std(data) #standard deviation of data
        prec = 1/sd #precision of the data
        res = array([exp(like.Normal(data,mu,prec)) for mu in arange(ll,ul,step)])  
        lik = res/max(res) # Likelihood function   
        print max(lik), min(lik)
    elif dist == 'exponential':
        res = [lamb**n*exp(-lamb*sum(data)) for lamb in arange(ll,ul,step)]
        lik = array(res)/max(array(res))
 
    elif dist == 'bernoulli':
        if ll<0 or ul>1:
            print "Parameter p of the bernoulli is out of range[0,1]"
        res = [exp(like.Bernoulli(data,p)) for p in arange(ll,ul,step)]
        lik = array(res)/max(array(res))
        
    elif dist == 'poisson':
        res = [exp(like.Poisson(data,lb)) for lb in arange(ll,ul,step)]
        lik = array(res)/max(array(res))
    
    elif dist == 'lognormal':
        sd = std(data) #standard deviation of data
        prec = 1/sd #precision of the data
        res = [exp(like.Lognormal(data,mu,prec)) for mu in arange(ll,ul,step)]
        lik = array(res)/max(array(res))    
    else:
        print 'Invalid distribution type. Valid distributions: normal, exponential, bernoulli and poisson'
    return lik


def Filt(cond, x, (ll, ul)):
    """
    filtering out Out-of-boundary thetas and phis. 
    for single output models.
    ul and ll are the pre-model boundaries of phi. 
    cond is a vector over which the conditional operations will be applied. 
    x is a vector or matrix of data. matrices are filtered line by line
    """
    #print cond.shape, x.shape, ll, ul
    cond = array(cond)
    cond = cond.ravel()
    if isinstance(x,tuple):
        l = len(x)
        x = array(x)
        x.shape = (l,x.size/float(l))
        #print 'shape of x is', x.shape
    else:
        #print 'shape of x is', x.shape
        pass
    try:
        f = compress(less(cond,ul) & greater(cond,ll),x, axis=1)
    except:
        f = compress(less(cond,ul) & greater(cond,ll),x)
    
        
    return f

def FiltM(cond,x,limits):
    """
    Multiple condition filtering.
    for multiple output models
    cond is an array of condition vectors
    limits is a list of tuples (ll,ul) with the length of cond
    """
    cond = array(cond)
    cnd = ones(cond.shape[1],int)
    for i,j in zip(cond,limits):
        ll = j[0]
        ul = j[1]
        #print cond.shape,cnd.shape,i.shape,ll,ul
        cnd = cnd & less(i,ul) & greater(i,ll)
    f = compress(cnd,x, axis=1)
    return f
        

def SIR(alpha,q2phi,limits,q2type,q1theta, phi,L, lik=[]):
    """
    Sampling Importance Resampling.

    :Parameters:
        - `alpha`: pooling weight;
        - `q2phi`: premodel of phi(tuple of vectors);
        - `limits`: limits for q2phi (list/tuple of tuples);
        - `q2type`: dist. type of q2phi (list of strings);
        - `q1theta`: premodel dists of thetas (tuple);
        - `phi`: model output (tuple of vectors);
        - `L`: size of the resample.
        - `lik`: list of likelihoods available
    """
##==On Uniform Priors we have to trim the density borders========================
##  The Density estimation with a gaussian kernel, extends beyond the limits of
##  an uniform distribution, due to this fact, we clip the ends of the kde
##  output in order to avoid artifacts.
##===============================================================================
    np = len(q1theta) # Number of parameters(theta) in the model
    no = len(phi) #Number of output variables

    q2pd =[]
    for i in xrange(no):
        (ll,ul) = limits[i] # limits of q2phi[i]
        if q2type[i] == 'uniform':
            q2pd.append(KDE(q2phi[i],(ll,ul)))
        else:
            q2pd.append(KDE(q2phi[i]))
    q2phi = q2pd
#---filtering out Out-of-boundary thetas and phis-------------------------------

    phi_filt=[]
    print "shape de q1theta[0]: ",q1theta[0].shape
    q1theta2 = array(q1theta) #Temporary copy to allow multiple filtering

    phi_filt = FiltM(phi,phi,limits) #filter Phis
    #print type(phi_filt)
    if not phi_filt.any():
        print "Due to bad specification of the prior distributions or of the model\nthe inference can't continue. please verify that your priors include at least\npart of the range of the output variables."
        return None
    #Remove thetas that generate out-of-bound phis for every phi
    q1theta_filt = FiltM(phi,q1theta2,limits)
#    print "shape de q1theta_filt (ln272): ",q1theta_filt.shape
    q1theta2 = q1theta_filt

    phi_filt = array(phi_filt)
# TODO: check to see if thetas or phis get empty due to bad priors!!!!
#-------------------------------------------------------------------------------

#---Calculate Kernel Density of the filtered phis-----------------------------------------------------------------------
    q1ed = []
    for i in xrange(no):
        (ll,ul) = limits[i] # limits of q2phi[i]
        if q2type[i] == 'uniform':
#            print sum(isinf(phi_filt))
            q1ed.append(KDE(phi_filt[i],(ll,ul)))
        else:
            q1ed.append(KDE(phi_filt[i]))
    q1est = q1ed
#-------------------------------------------------------------------------------

##==============================================================================
##Now, the two priors for Phi q2phi (derived from prior information and q1est 
##(generated by the model from the q1theta(priors on the inputs)), are pooled.
##The pooling is done by logarithmic pooling using alpha as a weighting factor.
##The higher the value of alpha the more wight is given to q1est.
##==============================================================================
#---Calculating the pooled prior of Phi-----------------------------------------
    qtilphi = []
    for i in xrange(no):
        qtilphi.append((array(q2phi[i]['y'])**(1-alpha))*(array(q1est[i]['y'])**alpha))
    qtilphi = array(qtilphi)
#-------------------------------------------------------------------------------
#---Calculating first term of the weigth expression-----------------------------
# TODO: Consider having a different alpha for each phi
    denslist=[]
    for i in xrange(no):
        #pairwise pooling of the phis and q2phis
        denslist.append((array(q2phi[i]['y'])/array(q1est[i]['y']))**(1-alpha))

    firstterm = denslist#product(denslist, axis=0)
#---Weights---------------------------------------------------------------------
        
    if not lik:
        w = firstterm #---- without likelihoods -----# 
    else:
        if len(lik)>1:
            prodlik = product(array(lik),axis=0)
        else:
            #only one likelihood function
            prodlik = lik[0]
#        w = firstterm*prodlik
        w = [i*prodlik for i in firstterm]
#-------------------------------------------------------------------------------
##========Link weights with each phi[i]=========================================
##  The weight vector (w) to be used in the resampling of the thetas is calculated
##  from operations on  densities. Consequently,its values are associated with
##  values on the support of Phi, not with the actual Phi[i] as output by the
##  model. Thus, its is necessary to recover the association between 
##  the Phi[i] (the outputs of each model run), and the weights
##  associated with them. For that, the support for phi is divided into 1024 bins
##  (the length of the weight vector), and the filtered Phi[i] are assigned to these bins
##  according to their value. This mapping is represented by the variable phi_bins
##  in which each element is the bin number of the correponding element in Phi.
##  A new weight vector(wi) is then created in which the elements of w are posi-
##  tioned according to the position of the Phi[i] to which it corresponds. That
##  is: w[i] = w[phi_bin[i]] repeated for each element i.
##==============================================================================
    
    bin_bound = []
    phi_bins = []
    wi = []
    for i in xrange(no):
        (ll,ul) = limits[i] #limits of phi
        step = (ul-ll)/1024.
        bin_bound.append(arange(ll,ul,step)) # Bin boundaries of the weight vector
        phi_bins.append(searchsorted(bin_bound[i], phi_filt[i])) # Return a vector of the bins for each phi
    g = lambda x:w[i][x-1]   # searchsorted returns 1 as the index for the first bin, not 0
    phi_bins = array(phi_bins)
    for i in xrange(no):
        wi.append(map(g,phi_bins[i]))
    wi = mean(array(wi),axis=0) #ATTENTION: Should this be averaged?

##========Resampling q1theta=====================================================
##  Here, the filtered q1theta are resampled according to the weight vector.  
##  L values are generated as indices to the weight vector wi(resamples) and used to resample
##  the parameters.
##===============================================================================

    # A given value is going to be resampled if random() < wi
    # A column of q1theta_filt is extracted for each value in resamples
    q = [0]*L
    wi = nan_to_num(array(wi))
    print sum(wi)
    if sum(wi) == 0:
        sys.exit('Resampling weights are all zero, please check your model or data.')
    j = 0
    while j < L: # Extract L samples from q1theta_filt
        i=randint(0,wi.size)# Random position of wi and q1theta_filt
        if random()<= wi[i]: 
            #print i, q1theta_filt.shape
            q[j]=q1theta_filt[:,i]# retain the sample according with resampling prob.
            j+=1
    # q is a list of arrays which is converted to an array and then transposed.
    #print "shape de q",len(q),q[0].shape
    qtiltheta = transpose(array(q)) 
    #print qtiltheta.shape
    return (w, qtiltheta, qtilphi, q1est)



# TODO: Implement calculation of Bayes factors!
#-------------------------------------------------------------------------------
##==MAIN========================================================================
#-------------------------------------------------------------------------------

def plotRaHist(arr):
    '''
    Plots a record array
    as a panel of histograms
    '''
    nv = len(arr.dtype.names)
    fs = (ceil(sqrt(nv)),floor(sqrt(nv))+1) #figure size
    P.figure()
    for i,n in enumerate(arr.dtype.names):
        P.subplot(nv/2+1,2,i+1)
        P.hist(arr[n],bins=50, normed=1, label=n)
        P.legend()

def main():
    """
    testing function
    """
    start = time()
    k = 20000 # Number of model runs
    L = 2000
    ll = 6
    ul = 9
    #data = [7,8,7,8,7,8,7]
    data = normal(7.5,1,400)
    lik = [] #initialize list of likelihoods
    lik.append(Likeli(data,'normal',(ll,ul)))
    
    q2phi = lhs.lhs(stats.uniform, (ll, ul), k)
    
    (phi, q1theta) = Run(k) # Runs the model
    print len(q1theta)
#---Restricting the range of phi------------------------------------------------
    
    (w, post_theta, qtilphi, q1est) = SIR(0.5,[q2phi],[(ll,ul)], ['uniform'],q1theta, [phi],L, lik)
    print "out of SIR"
    print post_theta.shape
#--generating the posterior of phi-------------------------------------------------------
    r = randint(0,len(post_theta[0]),L) #random index for the marginal posterior of r
    p = randint(0,len(post_theta[1]),L) #random index for the marginal posterior of p0
    post_phi = zeros(L,float) #initializing post_phi
    for i in xrange(L): #Monte Carlo with values of the posterior of Theta
        post_phi[i] = model(post_theta[0][r[i]],post_theta[1][p[i]])[-1]

    end = time()
    print end-start, ' seconds'
#---Plotting with matplotlib----------------------------------------------------------------------------
    P.figure(1)
    P.subplot(411)
    P.hist(post_theta[0],bins=50)
    P.ylabel(r'$\pi^{[r]}(\theta)$',fontsize=18)
    P.title('Posteriors and weight vector')
    P.subplot(412)
    P.hist(post_theta[1],bins=50)
    P.ylabel(r'$\pi^{[P_0]}(\theta)$',fontsize=18)
    P.subplot(413)
    P.hist(post_phi,bins=50)
    P.ylabel(r'$\pi^{[P]}(\phi)$',fontsize=18)
    ##plot(q1est['x'],qtilphi)
    ##ylabel(r'$P$', fontsize=12)
    P.subplot(414)
    P.plot(w)
    P.ylabel(r'$W_i$', fontsize=12)
    
    
    P.figure(2)
    P.subplot(411)
    P.hist(q1theta[0],bins=50)
    P.ylabel(r'$\theta r$',fontsize=18)
    P.title('Priors')
    P.subplot(412)
    P.hist(phi,bins=50)
    P.ylabel(r'$\phi$',fontsize=18)
    P.subplot(413)
    P.hist(q1theta[1],bins=50)
    P.ylabel(r'$\theta p_0$',fontsize=18)
    P.subplot(414)
    P.hist(q2phi,bins=50)
    P.ylabel(r'$q_2 \phi$',fontsize=18)
    P.show()

def main2():
    start = time()
    Me = Meld(K=10000,L=2000,model=model, ntheta=2,nphi=1,verbose=True)
    Me.setTheta(['r','p0'],[stats.uniform,stats.uniform],[(2,4),(0,5)])
    Me.setPhi(['p'],[stats.uniform],[(6,9)],[(6,9)])
    #Me.addData(normal(7.5,1,400),'normal',(6,9))
    #Me.run()
    Me.sir(data ={'p':[7.5]} )
    pt,pp = Me.getPosteriors()
    end = time()
    plotRaHist(pt)
    plotRaHist(pp)
    P.show()
    print end-start, ' seconds'
    
if __name__ == '__main__':
#    main()
    main2()
     

