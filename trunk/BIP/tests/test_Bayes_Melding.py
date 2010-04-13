import unittest
from BIP.Bayes import Melding
from numpy.testing import *
import scipy.stats as st

class TestFitModel(unittest.TestCase):
    def test___init__(self):
        K=50
        L=10
        ntheta = 2
        nphi = 3
        inits = [.999,0.001,0]
        tf = 25
        thetanames = ['beta','tau']
        phinames = ['S','I','R']
        verbose = False
        
        def model(*theta):
            try:#for standalone runs
                inits = y0
            except NameError:#for runing within a FitModel object
                inits = self.inits
            step = .1
            beta,tau = theta
            def sir(y,t):
                '''ODE model'''
                S,I,R = y
                return  [-beta*I*S, #dS/dt
                        beta*I*S - tau*I, #dI/dt
                        tau*I, #dR/dt
                        ]
            y = odeint(sir,inits,np.arange(t0,tf,step))
            return y
        self.model = model
        self.fit_model = Melding.FitModel(K, L, model, ntheta, nphi, inits, tf, thetanames, phinames, verbose)
        assert isinstance(self.fit_model,Melding.FitModel)

    def test_do_inference(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.do_inference(prior, data, predlen, method))
        assert True # TODO: implement your test here

    def test_init_priors(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.init_priors(prior))
        assert True # TODO: implement your test here

    def test_monitor(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.monitor())
        assert True # TODO: implement your test here

    def test_plot(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.plot())
        assert True # TODO: implement your test here

    def test_run(self):
        self.test_set_priors()
        wl = None
        nw = 1
        d = self.model(*[1.4,0.5])
        data = {'S':d[:,0],'I':d[:,1],'R':d[:,2]}
        method = 'SIR'
        monitor = False
        self.fit_model.run(wl, nw, data, method, monitor)
        assert self.fit_model.done_running 

    def test_plot_results(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, wl, nw, verbose)
        # self.assertEqual(expected, fit_model.plot_results(names))
        assert True # TODO: implement your test here

    def test_set_priors(self):
        tdists = [st.norm]*2
        tpars = [(1.4,.2),(0.5,.1)]
        tlims = [(0,3),(0,1)]
        pdists = [st.uniform]*3
        ppars = [(0,1)]*3
        plims = [(0,1)]*3
        self.fit_model.set_priors(tdists, tpars, tlims, pdists, ppars, plims)
        # self.assertEqual(expected, fit_model.set_priors(tdists, tpars, tlims, pdists, ppars, plims))
        assert True # TODO: implement your test here

class TestMeld(unittest.TestCase):
    def test___init__(self):
        self.meld = Meld(2000, 100, Melding.model, 2, 1, alpha=0.5, verbose=False, viz=False)
        assert isinstance(self.meld, Melding.Meld) 

    def test_abcRun(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.abcRun(fitfun, data, t, nopool, savetemp))
        assert False # TODO: implement your test here

    def test_add_salt(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.add_salt(dataset, band))
        assert False # TODO: implement your test here

    def test_basicfit(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.basicfit(s1, s2))
        assert False # TODO: implement your test here

    def test_filtM(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.filtM(cond, x, limits))
        assert False # TODO: implement your test here

    def test_getPosteriors(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.getPosteriors(t))
        assert False # TODO: implement your test here

    def test_imp_sample(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.imp_sample(n, data, w))
        assert False # TODO: implement your test here

    def test_logPooling(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.logPooling(phi))
        assert False # TODO: implement your test here

    def test_mcmc_run(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.mcmc_run(data, t, likvariance, burnin, nopool, method))
        assert False # TODO: implement your test here

    def test_mh(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.mh(n, t, data, likfun, variance, burnin))
        assert False # TODO: implement your test here

    def test_run(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.run(*args))
        assert False # TODO: implement your test here

    def test_runModel(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.runModel(savetemp, t, k))
        assert False # TODO: implement your test here

    def test_setPhi(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.setPhi(names, dists, pars, limits))
        assert False # TODO: implement your test here

    def test_setPhiFromData(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.setPhiFromData(names, data, limits))
        assert False # TODO: implement your test here

    def test_setTheta(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.setTheta(names, dists, pars))
        assert False # TODO: implement your test here

    def test_setThetaFromData(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.setThetaFromData(names, data, limits))
        assert False # TODO: implement your test here

    def test_sir(self):
        self.meld.sir(data ={'p':[7.5]} )
        pt,pp = self.meld.getPosteriors(1)
        # self.assertEqual(expected, meld.sir(data, t, variance, nopool, savetemp))
        assert_almost_equal(pp.P.mean(),7.5,1)

class TestEnumRun(unittest.TestCase):
    def test_enum_run(self):
        # self.assertEqual(expected, enumRun(model, theta, k))
        assert False # TODO: implement your test here

class TestModel(unittest.TestCase):
    def test_model(self):
        # self.assertEqual(expected, model(r, p0, n))
        assert False # TODO: implement your test here

class TestPlotRaHist(unittest.TestCase):
    def test_plot_ra_hist(self):
        # self.assertEqual(expected, plotRaHist(arr))
        assert False # TODO: implement your test here

class TestMain2(unittest.TestCase):
    def test_main2(self):
        # self.assertEqual(expected, main2())
        assert False # TODO: implement your test here

class TestMhTest(unittest.TestCase):
    def test_mh_test(self):
        # self.assertEqual(expected, mh_test())
        assert False # TODO: implement your test here

if __name__ == '__main__':
    unittest.main()