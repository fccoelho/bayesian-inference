from numpy.random import uniform, multinomial, exponential
#from numpy import arange, array, empty,zeros
import numpy as np
cimport numpy as np
import time
from random import random

DTYPE = np.double
ctypedef np.double_t DTYPE_t
ctypedef np.int_t INT_t

cdef extern from "math.h":
    double log(double)
cdef double clog(double x):
    return log(x*x)

cdef extern from "stdlib.h":
    double RAND_MAX
    double c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)
cdef crandom():
    cdef double rm = RAND_MAX
    return c_libc_random()/rm

cdef class Model(object):
    cdef object vn,pv, inits
    cdef np.ndarray tm,res,time,series, rates
    cdef unsigned int pvl,nvars,steps
    cdef object ts
    def __init__(self,vnames,rates,inits, tmat,propensity):
        '''
         * vnames: list of strings
         * rates: list of fixed rate parameters
         * inits: list of initial values of variables
         * propensity: list of lambda functions of the form: 
            lambda r,ini: some function of rates ans inits.
        '''
        self.vn = vnames
        self.rates = np.array(rates)
        self.inits = tuple(inits)
        self.tm = tmat
        self.pv = propensity#[compile(eq,'errmsg','eval') for eq in propensity]
        self.pvl = len(self.pv) #length of propensity vector
        self.nvars = len(self.inits) #number of variables
        self.time = np.zeros(1)
        self.series = np.zeros(1)
        self.steps = 0

    cpdef run(self, method='SSA', int tmax=10, int reps=1):
        cdef np.ndarray[np.int_t,ndim=3] res = np.zeros((tmax,self.nvars,reps),dtype=np.int)
        cdef np.ndarray[np.int_t] tvec = np.arange(tmax,dtype=np.int)
        self.res = res
        cdef int i, steps
        if method =='SSA':
            for i from 0 <= i<reps:
                steps = self.GSSA(tmax,i)
            #print steps,' steps'
        elif method == 'SSAct':
            pass
        self.time=tvec
       # self.series=self.res
        self.steps=steps
    
    def getStats(self):
        return self.time,self.res,self.steps
    
    cdef int GSSA(self, unsigned int tmax=50, unsigned int round=0) except *:
        '''
        Gillespie Direct algorithm
        '''
        cdef np.ndarray[np.int_t] ini = np.array(self.inits)
        cdef np.ndarray[np.double_t] r = self.rates
        pvi = self.pv
        cdef int l,steps,i,tim
        cdef double tc, tau, a0
        #cdef np.ndarray[INT_t] tvec
        cdef np.ndarray[np.double_t] pv 
        l=self.pvl
        pv = np.zeros(l,dtype=np.double)
        tm = self.tm
        #tvec = np.arange(tmax,dtype=int)
        tc = 0
        steps = 0
        self.res[0,:,round]= ini
        a0=1.
        for tim from 1<= tim <tmax:
            while tc < tim:
                for i from 0 <= i <l:
                    pv[i] = pvi[i](r,ini)
                #pv = abs(array([eq() for eq in pvi]))# #propensity vector
                a0 = a_sum(pv,l) #sum of all transition probabilities
                #print ini#,tim, pv, a0
                tau = (-1/a0)*clog(crandom())
                if pv.any():
                    event = np.random.multinomial(1,(pv/a0)) # event which will happen on this iteration
                    ini += tm[:,event.nonzero()[0][0]]
                #print tc, ini
                tc += tau
                steps +=1
                if a0 == 0: break
            self.res[tim,:,round] = ini
            if a0 == 0: break
#        tvec = tvec[:tim]
#        self.res = res[:tim,:,round]
        return steps

    def CR(self,pv):
        """
        Composition reaction algorithm
        """
        pass
        
cpdef double l1(np.ndarray r,np.ndarray ini):
    return r[0]*ini[0]*ini[1]
cpdef double l2(np.ndarray r,np.ndarray ini):
    return r[1]*ini[1]

cpdef main():
    vars = ['s','i','r']
    cdef np.ndarray[np.int_t] ini= np.array([500,1,0],dtype = np.int)
    cdef np.ndarray[np.float64_t] rates = np.array([.001,.1],dtype=np.float64)
    cdef np.ndarray[np.int_t,ndim=2] tm = np.array([[-1,0],[1,-1],[0,1]],dtype=np.int)
    prop = [l1,l2]
    M = Model(vnames = vars,rates = rates,inits=ini, tmat=tm,propensity=prop)
    t0=time.time()
    M.run(tmax=80,reps=1000)
    print 'total time: ',time.time()-t0
    return M.getStats()

    
cdef double a_sum(np.ndarray a, unsigned int len):
    cdef double s
    cdef int i
    s = 0
    for i from 0 <= i < len:
        s += a[i]
    return s



