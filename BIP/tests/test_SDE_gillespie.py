from __future__ import absolute_import
from __future__ import print_function
import nose
from nose import SkipTest
from nose.tools import assert_equal
from numpy import array,  zeros
from numpy.testing import assert_array_equal
from BIP.SDE.gillespie import Model

class TestDispatch:
    def test_dispatch(self):
        # self.assertEqual(expected, dispatch(model))
        raise SkipTest # TODO: implement your test here

def p1(r,ini): return r[0]*ini[0]*ini[1]
def p2(r,ini): return r[1]*ini[1]

class TestModel:
    def setUp(self):
        '''
        Create model for tests
        '''
        vnames = ['S','I','R']
        self.ini= [500,1,0]
        rates = [.001,.1]
        self.tm = array([[-1,0],[1,-1],[0,1]])
        #prop=[lambda r, ini:r[0]*ini[0]*ini[1],lambda r,ini:r[0]*ini[1]]
        self.M = Model(vnames = vnames,rates = rates,inits=self.ini, tmat=self.tm,propensity=[p1,p2])
    
    def test_shape_of_series(self):
        '''
        Tests if shape of series matrix is as expected
        '''
        model = self.M
        model.run(tmax=80,reps=10,viz=0,serial=1)
        t,series,steps, evts = model.getStats()
        assert_equal((10, 80, 3), series.shape) #(reps,t,vars)
    
    def test_events_match_series(self):
        '''
        Count number of events per time unit and make sure they correspond to change in time series 
        '''
        model = self.M
        model.run(tmax=80,reps=1,viz=0,serial=1)
        t,series,steps, evts = model.getStats()
        events = evts[0] #first replicate
        for n, ti in enumerate(t[1:]):
            print(ti)
            nev = 0
            ecv = zeros(3)
            for e, s in list(events.items()):
                evs = sum((s<=ti) & (s>t[n])) 
                nev += evs#number events of type e happening up until ti and after previous time-step
                print("n:", n,"t:", ti,  "nev:", evs, "evtype:", e,  series[0, n, :], series[0, n+1, :],  "evtimes:", s[:nev+1])
                if e ==0: #these events reduce by one the first variable
                    assert_equal(evs, series[0, n, 0]-series[0, n+1, 0])
                elif e == 1:
                    assert_equal(-evs, series[0, n, 2]-series[0, n+1, 2])
            ecv += abs(evs * self.tm[:, e]) #expected change in variables
        
        
    def test_model(self):
        # assert_equal(expected, model(r, p0, n))
        raise SkipTest # TODO: implement your test here

    def TestGSSA(self):
        # model = Model(vnames, rates, inits, tmat, propensity)
        # assert_equal(expected, model.GSSA())
        raise SkipTest # TODO: implement your test here

    def test_getStats(self):
        model = self.M
        model.run(tmax=80,reps=10,viz=0,serial=1)
        assert isinstance(model.getStats(), tuple)
        
    def test_run_serial(self):
        model = self.M
        model.run(tmax=80,reps=10,viz=0,serial=1)
    
    def test_run_parallel(self):
        model = self.M
        model.run(tmax=80,reps=10,viz=0,serial=0)
        

class TestP1:
    def test_p1(self):
        # self.assertEqual(expected, p1(r, ini))
        raise SkipTest # TODO: implement your test here

class TestP2:
    def test_p2(self):
        # self.assertEqual(expected, p2(r, ini))
        raise SkipTest # TODO: implement your test here

class TestMain:
    def test_main(self):
        # assert_equal(expected, main())
        raise SkipTest # TODO: implement your test here

if __name__ == '__main__':
    nose.run()
