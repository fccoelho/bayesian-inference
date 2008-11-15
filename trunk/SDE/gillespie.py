from numpy.random import uniform, multinomial, exponential,random
from numpy import arange, array, empty,zeros,log
#from math import log
import time
from multiprocessing import Pool



import psyco
psyco.full()

#global ini
#ini=[500,1,0]

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
        self.nvars = len(self.inits) #number of variables
        self.time=None
        self.series=None
        self.steps=0
    
    def getStats(self):
        return self.time,self.series,self.steps
    
    def run(self,method='SSA', tmax=10, reps=1):
        self.res = zeros((tmax,self.nvars,reps),dtype=float)
        tvec = arange(tmax, dtype=int)
        pool = Pool()    
        if method =='SSA':
            #r = pool.imap(self.GSSA, ((tmax,i) for i in xrange(reps)))
            #steps = r.get()
            for i in xrange(reps):
                steps = self.GSSA(tmax,i)
            print steps,' steps'
        elif method == 'SSAct':
            pass
        self.time=tvec
        self.series=self.res
        self.steps=steps
    def GSSA(self, tmax=50,round=0):
        '''
        Gillespie Direct algorithm
        '''
        ini = self.inits
        r = self.rates
        pvi = self.pv
        l=self.pvl
        pv = zeros(l,dtype=float)
        tm = self.tm
        
        tc = 0
        steps = 0
        self.res[0,:,round]= ini
        a0=1
        for tim in xrange(1,tmax):
            while tc < tim:
                for i in xrange(l):
                    pv[i] = pvi[i](r,ini)
                #pv = abs(array([eq() for eq in pvi]))# #propensity vector
                a0 = pv.sum() #sum of all transition probabilities
#                print tim, pv, a0
                tau = (-1/a0)*log(random())
                event = multinomial(1,pv/a0) # event which will happen on this iteration
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
        

def main():
    vars = ['s','i','r']
    ini= [500,1,0]
    rates = [.001,.1]
    tm = array([[-1,0],[1,-1],[0,1]])
    
    prop=[lambda r, ini:r[0]*ini[0]*ini[1],lambda r,ini:r[0]*ini[1]]
    M = Model(vnames = vars,rates = rates,inits=ini, tmat=tm,propensity=prop)
    t0=time.time()
    M.run(tmax=80,reps=1000)
    print 'total time: ',time.time()-t0
    #print res

#    from pylab import plot , show, legend
#    plot(t,res,'-.')
#    legend(M.vn,loc=0)
#    show()
    

if __name__=="__main__":
    import cProfile
    cProfile.run('main()',sort=1)
    main()
