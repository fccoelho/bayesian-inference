# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
# Name:        gillespie.py
# Project:	Bayesian-Inference
# Purpose:     
#
# Author:      Flávio Codeço Coelho<fccoelho@gmail.com>
#
# Created:     2008-11-26
# Copyright:   (c) 2008 by the Author
# Licence:     GPL
#-----------------------------------------------------------------------------
__docformat__ = "restructuredtext en"
from numpy.random import uniform, multinomial, exponential,random
from numpy import arange, array, empty,zeros,log
#from math import log
import time
from multiprocessing import Pool



import psyco
psyco.full()

#global ini
#ini=[500,1,0]
def dispatch(model):
    '''this function is necessary for paralelization'''
    return model.GSSA()

class Model:
    def __init__(self,vnames,rates,inits, tmat,propensity):
        '''
         * vnames: list of strings
         * rates: list of fixed rate parameters
         * inits: list of initial values of variables
         * propensity: list of lambda functions of the form: 
            lambda r,ini: some function of rates ans inits.
        '''
        self.vn = vnames
        self.rates = rates
        self.inits = inits
        self.tm = tmat
        self.pv = propensity#[compile(eq,'errmsg','eval') for eq in propensity]
        self.pvl = len(self.pv) #length of propensity vector
        self.pv0 = zeros(self.pvl,dtype=float)
        self.nvars = len(self.inits) #number of variables
        self.time=None
        self.tmax = None
        self.series=None
        self.steps=0
    
    def getStats(self):
        return self.time,self.series,self.steps
    
    def run(self,method='SSA', tmax=10, reps=1):
        self.tmax = tmax
        #self.res = zeros((tmax,self.nvars,reps),dtype=float)
        self.res = zeros((tmax,self.nvars),dtype=float)
        tvec = arange(tmax, dtype=int)
        pool = Pool()    
        if method =='SSA':
            # Parallel version
            self.res = array(pool.map(dispatch,[(self) for i in xrange(reps)],chunksize=10))
            # Serial
            #self.res = array([self.GSSA() for i in xrange(reps)])
            
        elif method == 'SSAct':
            pass
        pool.close()
        pool.join()
        self.time=tvec
        self.series=self.res
        #self.steps=steps
    
    def GSSA(self):
        '''
        Gillespie Direct algorithm
        '''
        tmax = self.tmax
        ini = self.inits
        r = self.rates
        pvi = self.pv #propensity functions
        tm = self.tm
        pv = self.pv0 #actual propensity values for each time step
        tc = 0
        steps = 0
        self.res[0,:]= ini
        for tim in xrange(1,tmax):
            while tc < tim:
                i=0
                a0=0
                for p in pvi:
                    pv[i] = p(r,ini)
                    a0+=pv[i]
                    i+=1
                tau = (-1/a0)*log(random())
                event = multinomial(1,pv/a0) # event which will happen on this iteration
                ini += tm[:,event.nonzero()[0][0]]
                #print tc, ini
                tc += tau
                steps +=1
                if a0 == 0: break
            self.res[tim,:] = ini
            if a0 == 0: break
        return self.res


def p1(r,ini): return r[0]*ini[0]*ini[1]
def p2(r,ini): return r[0]*ini[1]

def main():
    vars = ['s','i','r']
    ini= [500,1,0]
    rates = [.001,.1]
    tm = array([[-1,0],[1,-1],[0,1]])
    #prop=[lambda r, ini:r[0]*ini[0]*ini[1],lambda r,ini:r[0]*ini[1]]
    M = Model(vnames = vars,rates = rates,inits=ini, tmat=tm,propensity=[p1,p2])
    t0=time.time()
    M.run(tmax=80,reps=1000)
    print 'total time: ',time.time()-t0
    #print res

#    from pylab import plot , show, legend
#    plot(t,res,'-.')
#    legend(M.vn,loc=0)
#    show()
    

if __name__=="__main__":
    #import cProfile
    #cProfile.run('main()',sort=1,filename='gillespie.profile')
    main()
