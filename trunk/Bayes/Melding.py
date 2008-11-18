# -*- coding:latin-1 -*-
#-----------------------------------------------------------------------------
# Name:        Melding.py
# Purpose:     Bayesian melding
#
# Author:      Flávio Codeço Coelho
#
# Created:     2003/08/10
# RCS-ID:      $Id: Melding.py $
# Copyright:   (c) 2003
# Licence:     GPL
# New field:   Whatever
#-----------------------------------------------------------------------------
##import psyco
##psyco.full()
import scipy.stats.kde as kde
#from random import *
#from cmath import *
import pylab as P
from numpy import *
from numpy.random import *
import lhs
#import random
import like
import sys



def model(r, p0, n=1):
    """
    Model (r,p0, n=1)
    Simulates the Population dynamic Model (PDM) Pt = rP0
    for n time steps.
    P0 is the initial population size. 
    Example model for testing purposes.
    """
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
#---q1theta---------------------------------------------------------------------
    #r = genprior('uniform',(2,4), k) # Uniform [2,4]
    #r = lhs.lhs(['r'],['Uniform'],[[2,4]],k)[0]
    #p0 = genprior('uniform',(0,5), k)# Uniform[0,5]
    #p0 = lhs.lhs(['p0'],['Uniform'],[[0,5]],k)[0]
    r,p0 = lhs.lhs(['r', 'p0'],['uniform','uniform'],[[2,4],[0,5]],k,noCorrRestr=False)
    q1theta = (r, p0)
#-------------------------------------------------------------------------------
    phi=zeros(k, float)
    print r.shape, p0.shape
    for i in xrange(k):
        phi[i] = model(r[i], p0[i])[-1] # Sets phi[i] to the last point of the simulation
        
    
    return phi, q1theta

def KDE(x, (ll, ul)=('','')):
    """
    KDE(x)
    performs a kernel density estimate using the R function density
    if (ll,ul) enforce limits for the distribution.
    Returns a dictionary.
    """
    #r.assign("x", x)
    
    if not ll == '':
        rn=arange(ll,ul,(ul-ll)/512.)
        #print x.shape,rn.shape
        est = kde.gaussian_kde(x.ravel()).evaluate(rn)
        #r.assign("ll", ll)
        #r.assign("ul", ul)
        #est = r('density(x,from=ll, to=ul)') #trims the density borders
    else:
        ll = min(x)
        ul = max(x)
        rn=arange(ll,ul,(ul-ll)/512.)
        est = kde.gaussian_kde(x).evaluate(rn)
        #est = r('density(x)')
        print 'No - KDE'
        
    
    
    return {'y':est,'x':rn}


def Likeli(data, dist, limits,**kwargs):
    """
    Generates the likelihood function of data given dist.
    limits is a tuple setting the interval of the parameter space that will
    be used as the support for the Likelihood function.
    returns a vector (512 elements).
    """
    n = len(data) # Number of data points
    data = array(data)
    (ll,ul) = limits #limits for the parameter space
    step = (ul-ll)/512.
    
    if dist == 'normal': # In this case, L is a function of the mean.
        sd = std(data) #standard deviation of data
        prec = 1/sd #precision of the data
        res = [exp(like.Normal(data,mu,prec)) for mu in arange(ll,ul,step)] #empty list of results
        #for mu in arange(ll,ul,step):
            #print res#mu, sd
            #res.append(exp(-0.5*sum(((data-mu)/sd)**2))) # Makes a list of the likelihood of every datum. Removed this constant term from formula because of dependence on n :1./((sd*sqrt(2*pi))**n)
        res = array(res)    
        lik = res/max(res) # Likelihood function   
        print max(lik), min(lik)
    elif dist == 'exponential':
        res = [lamb**n*exp(-lamb*sum(data)) for lamb in arange(ll,ul,step)]
##        for lamb in arange(ll,ul,step):
##            res.append(lamb**n*exp(-lamb*sum(data)))
        lik = array(res)/max(array(res))
 
    elif dist == 'bernoulli':
        if ll<0 or ul>1:
            print "Parameter p of the bernoulli is out of range[0,1]"
        res = [exp(like.Bernoulli(data,p)) for p in arange(ll,ul,step)]
        #for p in arange(ll,ul,step):
            #res.append(p**sum(data)*(1-p)**(n-sum(data)))
            #res.append(exp(flib.bernoulli(data,p)))
        lik = array(res)/max(array(res))
        
    elif dist == 'poisson':
        res = [exp(like.Poisson(data,lb)) for lb in arange(ll,ul,step)]
##        for lb in arange(ll,ul,step):
##            res.append(exp(-n*lb)*(lb*sum(data))/product(factorial(data)))
        lik = array(res)/max(array(res))
    
    elif dist == 'lognormal':
        sd = std(data) #standard deviation of data
        prec = 1/sd #precision of the data
        res = [exp(like.Lognormal(data,mu,prec)) for mu in arange(ll,ul,step)]
        lik = array(res)/max(array(res))    
    else:
        print 'Invalid distribution type. Valid distributions: normal, exponential, bernoulli and poisson'
    
    return lik
    
def factorial(numb):
    """
    calculates the factorial of a number, or a sequence of numbers.
    """
    l = ["<type 'array'>","<type 'list'>","<type 'tuple'>"]
    if not str(type(numb)) in l:
        n = int(numb)
        if n == 0:
            return 1
        else:
            return n * factorial(n-1)
    else:
        res=[]
        for i in numb:
            res.append(factorial(i))
        return array(res)


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
    alpha = pooling weight; 
    q2phi = premodel of phi(vector or a tuple of vectors); 
    limits = limits for q2phi (tuple or list/tuple of tuples); 
    q2type = dist. type of q2phi (string or list of strings); 
    q1theta = premodel dists of thetas (tuple); 
    phi = model output (vector or tuple of vectors); 
    L = size of the resample.
    Lik = list of likelihoods available
    """
    
    
##==On Uniform Priors we have to trim the density borders========================
##  The Density estimation with a gaussian kernel, extends beyond the limits of
##  an uniform distribution, due to this fact, we clip the ends of the kde
##  output in order the avoid artifacts.
##===============================================================================
    np = len(q1theta) # Number of parameters(theta) in the model
#---for multicompartimental models-----------------------------------------------
    multi = ["<type 'list'>","<type 'tuple'>"] #possible types of data structures of q2phi and phi

    if not str(type(phi)) in multi:
        (ll,ul) = limits # limits of q2phi for single compartment models
    no = None
    if str(type(q2phi)) in multi:
        no = len(phi) #Number of output variables
        q2phi = array(q2phi)
        q2pd =[]
        for i in xrange(no):
            (ll,ul) = limits[i] # limits of q2phi[i]
            if q2type[i] == 'uniform':
                q2pd.append(KDE(q2phi[i],(ll,ul)))
            else:
                q2pd.append(KDE(q2phi[i]))
        q2phi = q2pd
#---for single compartiment models----------------------------------------------------------------------------   
    else:
        if q2type == 'uniform':
            q2phi = KDE(q2phi, (ll,ul)) #calculating the kernel density
            print 'Ok'
        else:
            q2phi = KDE(q2phi)
            print 'No - SIR1'
#-------------------------------------------------------------------------------
        
#---filtering out Out-of-boundary thetas and phis-------------------------------
    if str(type(phi)) in multi: #Multicompartimental models
        #phi = array(phi)# Convert list/tuple of vectors in array, where each vector becomes a line of the array.
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
        print "shape de q1theta_filt (ln272): ",q1theta_filt.shape
        q1theta2 = q1theta_filt
        
        ## for i in xrange(no):
##             (ll,ul) = limits[i] # limits of q2phi[i]
##             print 'no = %s'%no
##             phi_filt.append(Filt(phi[i],phi[i],(ll,ul))) #filter Phis
##             if not phi_filt[i].any():
##                 print "Due to bad specification of the prior distributions or of the model\nthe inference can't continue. please verify that your priors include at least\npart of the range of the output variables."
##                 return None
                
##             q1theta_filt = Filt(phi[i],q1theta2,(ll,ul)) #Remove thetas that generate out-of-bound phis for every phi
##             q1theta2 = q1theta_filt
            
        phi_filt = array(phi_filt)
# TODO: check to see if thetas or phis go to empty due to bad priors!!!!
    else: #Single compartment
        phi_filt = Filt(phi,phi,(ll,ul)) #remove out-of-bound phis
        q1theta_filt = Filt(phi,q1theta,(ll,ul)) #remove thetas that generate out-of-bound phis
#-------------------------------------------------------------------------------

#---Calculate Kernel Density of the filtered phis----------------------------------------------------------------------------
    if no: #multicompartimental
        q1ed = []
        for i in xrange(no):
            (ll,ul) = limits[i] # limits of q2phi[i]
            if q2type[i] == 'uniform':
                q1ed.append(KDE(phi_filt[i],(ll,ul)))
            else: 
                q1ed.append(KDE(phi_filt[i]))
        q1est = q1ed
    else: #Single compartment
        if q2type == 'uniform':
            q1est = KDE(phi_filt,(ll,ul)) # Density of of model outputs restricted to the range determined by the Priors, so that we can pool q1est and q2phi.
            print 'Ok'
        else:
            q1est = KDE(phi_filt)
            print 'No - SIR2'
#-------------------------------------------------------------------------------

##==============================================================================
##Now, the two priors for Phi q2phi (derived from prior information and q1est 
##(generated by the model from the q1theta(priors on the inputs)), are pooled.
##The pooling is done by logarithmic pooling using alpha as a weighting factor.
##The higher the value of alpha the more wight is given to q1est.
##==============================================================================
#---Calculating the pooled prior of Phi------------------------------------------
    if no: #multicompartimental
        qtilphi = []
        for i in xrange(no):
            qtilphi.append((array(q2phi[i]['y'])**(1-alpha))*(array(q1est[i]['y'])**alpha))
        qtilphi = array(qtilphi)
    else:  #Single compartment
        qtilphi = (array(q2phi['y'])**(1-alpha))*(array(q1est['y'])**alpha) # Pooled prior of Phi
#-------------------------------------------------------------------------------
    
    
#---Calculating first term of the weigth expression-----------------------------
# TODO: Consider having a different alpha for each phi
    if no:#multicompartimental
        denslist=[]
        for i in xrange(no):
            #pairwise pooling of the phis and q2phis
            denslist.append((array(q2phi[i]['y'])/array(q1est[i]['y']))**(1-alpha)) 
            
        firstterm = product(denslist, axis=0)
    else:#Single compartment
        firstterm = (array(q2phi['y'])/array(q1est['y']))**(1-alpha)
        
        
#---Weights---------------------------------------------------------------------
        
    if not lik:
        w = firstterm #---- without likelihoods -----# 
    else:
        if len(lik)>1:
            prodlik = product(array(lik),axis=0)
        else:
            #only one likelihood function
            prodlik = lik[0]
        w = firstterm*prodlik
         
#-------------------------------------------------------------------------------
##========Link weights with each phi[i]=========================================
##  The weight vector (w) to be used in the resampling of the thetas is calculated
##  from operations on  densities. Consequently,its values are associated with
##  values on the support of Phi, not with the actual Phi[i] as output by the
##  model. Thus, its is necessary to recover the asso-
##  ciation between the Phi[i] (the outputs of each model run), and the weights
##  associated with them. For that, the support for phi is divided into 512 bins
##  (the length of the weight vector), and the filtered Phi[i] are assigned to these bins
##  according to their value. This mapping is represented by the variable phi_bins
##  in which each element is the bin number of the correponding element in Phi.
##  A new weight vector(wi) is then created in which the elements of w are posi-
##  tioned according to the position of the Phi[i] to which it corresponds. That
##  is: w[i] = w[phi_bin[i]] repeated for each element i.
##==============================================================================
    
    if no:#multicompartimental
        bin_bound = []
        phi_bins = []
        wi = []
        for i in xrange(no):
            (ll,ul) = limits[i] #limits of phi
            step = (ul-ll)/512.     
            bin_bound.append(arange(ll,ul,step)) # Bin boundaries of the weight vector
            phi_bins.append(searchsorted(bin_bound[i], phi_filt[i])) # Return a vector of the bins for each phi
        g = lambda x:w[x-1]   # searchsorted returns 1 as the index for the first bin, not 0
        phi_bins = array(phi_bins)
        for i in xrange(no):
            wi.append(map(g,phi_bins[i]))
        wi = mean(array(wi),axis=0) #ATTENTION: Should this be averaged?
    else: #single compartment
        (ll,ul) = limits
        step = (ul-ll)/512.  
        bin_bound = arange(ll,ul,step) # Bin boundaries of the weight vector
        phi_bins = searchsorted(bin_bound, phi_filt) # Return a vector of the bins for each phi
        g = lambda x:w[x-1]   # searchsorted returns 1 as the index for the first bin, not 0
        wi = map(g,phi_bins.ravel())
#-------------------------------------------------------------------------------
   
        
#---creating a biased die based on probabilities of w---------------------------
    #die = cumsum(wi)#-Cumulative sum of resampling probabilities
    #roll = uniform(die[0],die[-1],L)
   
    
##========Resampling q1theta=====================================================
##  Here, the filtered q1theta are resampled according to the weight vector.  
##  L values are generated as indices to the weight vector wi(resamples) and used to resample
##  the parameters.
##===============================================================================
    #sampled_is = searchsorted(die, roll)
    #qtiltheta = transpose(array(map(h,sampled_is)))
    #qtiltheta=zeros((np,L), Float) # Initialize the qtiltheta matrix
    #resamples = randint(0,len(wi),(L,))# Random order of resamples

    # A given value is going to be resampled if random() < wi
    # A column of q1theta_filt is extracted for each value in resamples
    q = [0]*L
    wi = array(wi)
    if max(wi) ==0:
        sys.exit('Resampling weights are all zero, please check your model or data.')
    j = 0
    while j < L: # Extract L samples from q1theta_filt
        i=randint(0,wi.size)# Random position of wi and q1theta_filt
        if random()<= wi[i]: 
            #print i, q1theta_filt.shape
            q[j]=q1theta_filt[:,i]# retain the sample according with resampling prob.
            j+=1
    # q is a list of arrays which is converted to an array and then transposed.
    print "shape de q",len(q),q[0].shape
    qtiltheta = transpose(array(q)) 
    print qtiltheta.shape
    return (w, qtiltheta, qtilphi, q1est)

def plotmat(x, tit='title', b=50):
    """
    This funtion implements a simple 50 bin, normalized histogram using the matplotlib module.
    """

    P.hist(x,bins=b,normed=1)
    P.ylabel(tit, fontsize=18)

def genprior(type, params, shape=[]):
    """
    genprior(type, params, shape)
    The function returns a vector or a matrix containinin a sample of the specified distribution with size given by shape.
    """
    seed()
    distlist=['uniform', 'normal', 'exponential', 'beta', 'gamma', 'chi-square', 'F', 'binomial', 'neg-binomial', 'poisson', 'multinomial']
    if type == 'uniform':
        prior = uniform(params[0], params[1], shape)
    elif type == 'normal':
        prior = normal(params[0], params[1], shape)
    elif type == 'exponential':
        prior = exponential(params, shape)
    elif type == 'beta':
        prior = beta(params[0], params[1], shape)
    elif type == 'gamma':
        prior = gamma(params[0], params[1], shape)
    elif type == 'chi-square':
        prior = chi_square(params, shape)
    elif type == 'F':
        prior = F(params[0], params[1], shape)
    elif type == 'binomial':
        prior = binomial(params[0], params[1], shape)
    elif type == 'neg-binomial':
        prior = negative_binomial(params[0], params[1], shape)
    elif type == 'poisson':
        prior = poisson(params, shape)
    elif type == 'multinomial':
        prior = multinomial(params)
    else:
        print 'Invalid distribution type.'
    
    return prior


# TODO: Implement calculation of Bayes factors!
#-------------------------------------------------------------------------------
##==MAIN========================================================================
#-------------------------------------------------------------------------------


def main():
    """
    testing function
    """
    k = 20000 # Number of model runs
    L = 2000
    ll = 6
    ul = 9
    #data = [7,8,7,8,7,8,7]
    data = normal(7.5,1,400)
    lik = [] #initialize list of likelihoods
    lik.append(Likeli(data,'normal',(ll,ul)))
    #q2phi = genprior('uniform', (ll,ul), k) # uniform[6,9] - pre-model distribution of Phi
    q2phi = lhs.lhs(['p'],['uniform'],[[ll,ul]],k)[0]
    
    (phi, q1theta) = Run(k) # Runs the model
    print len(q1theta)
#---Restricting the range of phi------------------------------------------------
    
    (w, post_theta, qtilphi, q1est) = SIR(0.5,q2phi,(ll,ul), 'uniform',q1theta, phi,L, lik)
    print "out of SIR"
    print post_theta.shape
#--generating the posterior of phi-------------------------------------------------------
    r = uniform(0,len(post_theta[0]),L) #random index for the marginal posterior of r
    p = uniform(0,len(post_theta[1]),L) #random index for the marginal posterior of p0
    post_phi = zeros(L,float) #initializing post_phi
    for i in xrange(L): #Monte Carlo with values of the posterior of Theta
        post_phi[i] = model(post_theta[0][int(r[i])],post_theta[1][int(p[i])])[-1]

#---Plotting with matplotlib----------------------------------------------------------------------------
    P.figure(1)
    P.subplot(411)
    plotmat(post_theta[0], tit=r'$\pi^{[r]}(\theta)$')
    P.title('Posteriors and weight vector')
    P.subplot(412)
    plotmat(post_theta[1], tit=r'$\pi^{[P_0]}(\theta)$')
    P.subplot(413)
    plotmat(post_phi, tit=r'$\pi^{[P]}(\phi)$')
    ##plot(q1est['x'],qtilphi)
    ##ylabel(r'$P$', fontsize=12)
    P.subplot(414)
    P.plot(w)
    P.ylabel(r'$W_i$', fontsize=12)
    
    
    P.figure(2)
    P.subplot(411)
    plotmat(q1theta[0], tit=r'$\theta r$')
    P.title('Priors')
    P.subplot(412)
    plotmat(phi, tit=r'$\phi$')
    P.subplot(413)
    plotmat(q1theta[1], tit=r'$\theta p_0$')
    P.subplot(414)
    plotmat(q2phi, tit=r'$q_2 \phi$')
    P.show()
    
if __name__ == '__main__':
    from time import clock
    start = clock()
    main()  
    end = clock()
    print end-start, ' seconds'  
