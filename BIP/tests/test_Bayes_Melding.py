from __future__ import absolute_import
from __future__ import print_function
from BIP.Bayes import Melding
from numpy.testing import *
import scipy.stats as st
import numpy as np
from nose.tools import assert_equal, set_trace
from nose import SkipTest
from scipy.integrate import odeint
import unittest
import os, glob


K = 50
inits = [.999, 0.001, 0]
tf = 25
thetanames = ['beta', 'tau']
phinames = ['S', 'I', 'R']
verbose = 2


def model(*theta):
    step = .1
    beta, tau = theta

    def sir(y, t):
        '''ODE model'''
        S, I, R = y
        return [-beta * I * S,  # dS/dt
                beta * I * S - tau * I,  #dI/dt
                tau * I,  #dR/dt
                ]

    y = odeint(sir, inits, np.arange(0, tf, step))
    return y


fit_model = Melding.FitModel(K, model, inits, tf, thetanames, phinames, verbose=verbose)


class TestFitModel(unittest.TestCase):
    def setUp(self):
        nt = len(thetanames)
        np = len(phinames)
        fit_model.set_priors(tdists=nt * [st.uniform], tpars=[(0, 4), (0, 2)], tlims=[(0, 4), (0, 2)],
                             pdists=np * [st.uniform], ppars=np * [(0, 1)], plims=np * [(0, 1)])

    def tearDown(self):
        files = glob.glob('*.pickle') + glob.glob('*.sqlite')
        for f in files:
            os.unlink(f)

    def test___init__(self):
        assert isinstance(fit_model, Melding.FitModel)

    def test_check_priors(self):
        fit_model.check_priors()

    def test_do_inference(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.do_inference(prior, data, predlen, method))
        assert True  # TODO: implement your test here

    def test_init_priors(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.init_priors(prior))
        assert True  # TODO: implement your test here

    def test_monitor(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.monitor())
        assert True  # TODO: implement your test here

    def test_plot(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.plot())
        assert True  # TODO: implement your test here

    def test_run(self):
        self.test_set_priors()
        wl = None
        nw = 1
        d = model(*[1.4, 0.5])
        data = {'S': d[:, 0], 'I': d[:, 1], 'R': d[:, 2]}
        likvar = 1e-1
        method = 'SIR'
        fit_model.run(data, method, likvar, monitor=False)
        assert fit_model.done_running

    def test_plot_results(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, wl, nw, verbose)
        # self.assertEqual(expected, fit_model.plot_results(names))
        assert True  # TODO: implement your test here

    def test_set_priors(self):
        tdists = [st.norm] * 2
        tpars = [(1.4, .2), (0.5, .1)]
        tlims = [(0, 3), (0, 1)]
        pdists = [st.uniform] * 3
        ppars = [(0, 1)] * 3
        plims = [(0, 1)] * 3
        fit_model.set_priors(tdists, tpars, tlims, pdists, ppars, plims)
        # self.assertEqual(expected, fit_model.set_priors(tdists, tpars, tlims, pdists, ppars, plims))
        assert True  # TODO: implement your test here

    def TestAIC_from_RSS(self):
        # fit_model = FitModel(K, model, inits, tf, thetanames, phinames, wl, nw, verbose, burnin)
        # assert_equal(expected, fit_model.AIC_from_RSS())
        raise SkipTest  # TODO: implement your test here

    def test_optimize(self):
        # fit_model = FitModel(K, model, inits, tf, thetanames, phinames, wl, nw, verbose, burnin)
        # assert_equal(expected, fit_model.optimize(data, p0, optimizer, tol, verbose, plot))
        raise SkipTest  # TODO: implement your test here

    def test_prior_sample(self):
        # fit_model = FitModel(K, model, inits, tf, thetanames, phinames, wl, nw, verbose, burnin)
        # assert_equal(expected, fit_model.prior_sample())
        raise SkipTest  # TODO: implement your test here


class TestMeld:
    def setUp(self):
        self.meld = Melding.Meld(3000, 1000, Melding.model, 2, 1, alpha=0.5, verbose=0, viz=False)
        self.meld.setTheta(['r', 'p0'], [st.uniform, st.uniform], [(2, 4), (0, 5)], [(0, 10), (0, 10)])
        self.meld.setPhi(['p'], [st.uniform], [(6, 9)], [(0, 10)])

    def test___init__(self):
        assert isinstance(self.meld, Melding.Meld)

    def test_abcRun(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.abcRun(fitfun, data, t, nopool, savetemp))
        raise SkipTest  # TODO: implement your test here

    def test_add_salt(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.add_salt(dataset, band))
        raise SkipTest  # TODO: implement your test here

    def test_basicfit(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.basicfit(s1, s2))
        raise SkipTest  # TODO: implement your test here

    def test_filtM(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.filtM(cond, x, limits))
        raise SkipTest  # TODO: implement your test here

    def test_getPosteriors(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.getPosteriors(t))
        raise SkipTest  # TODO: implement your test here

    def test_imp_sample(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.imp_sample(n, data, w))
        raise SkipTest  # TODO: implement your test here

    def test_logPooling(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.logPooling(phi))
        raise SkipTest  # TODO: implement your test here

    def test_mcmc_run(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.mcmc_run(data, t, likvariance, burnin, nopool, method))
        raise SkipTest  # TODO: implement your test here

    def test_mh(self):
        self.meld.setTheta(['r', 'p0'], [st.uniform, st.uniform], [(2, 4), (0, 5)], [(0, 10), (0, 10)])
        self.meld.setPhi(['p'], [st.uniform], [(6, 9)], [(0, 10)])
        self.meld.mcmc_run(data={'p': np.array([7.5])}, burnin=500)
        pt, pp = self.meld.getPosteriors(1)
        assert_almost_equal(pp.p.mean(), 7.5, 1)

    def test_run(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.run(*args))
        raise SkipTest  # TODO: implement your test here

    def test_runModel(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.runModel(savetemp, t, k))
        raise SkipTest  # TODO: implement your test here

    def test_setPhi(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.setPhi(names, dists, pars, limits))
        raise SkipTest  # TODO: implement your test here

    def test_setPhiFromData(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.setPhiFromData(names, data, limits))
        raise SkipTest  # TODO: implement your test here

    def test_setTheta(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.setTheta(names, dists, pars))
        raise SkipTest  # TODO: implement your test here

    def test_setThetaFromData(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.setThetaFromData(names, data, limits))
        raise SkipTest  # TODO: implement your test here

    def test_sir(self):
        self.meld.sir(data={'p': [7.5]}, pool=True)
        # self.meld.run()
        pt, pp = self.meld.getPosteriors(1)
        # self.assertEqual(expected, meld.sir(data, t, variance, nopool, savetemp))
        assert_almost_equal(pp.p.mean(), 7.5, 1)

    def test_current_plot(self):
        d = np.sin(np.arange(100))
        series = np.recarray((1, 100), formats=['f8'], names=['sine'])
        series.sine[0, :] = d
        data = {'sine': d}
        fit_model.Me.current_plot(series, data, 0)
        fit_model.Me.every_run_plot.close_plot()


class TestEnumRun:
    def test_enum_run(self):
        # self.assertEqual(expected, enumRun(model, theta, k))
        raise SkipTest  # TODO: implement your test here


class TestMode:
    def test_model(self):
        # self.assertEqual(expected, model(r, p0, n))
        raise SkipTest  # TODO: implement your test here


class TestPlotRaHist:
    def test_plot_ra_hist(self):
        # self.assertEqual(expected, plotRaHist(arr))
        raise SkipTest  # TODO: implement your test here


class TestMain2:
    def test_main2(self):
        # self.assertEqual(expected, main2())
        raise SkipTest  # TODO: implement your test here


class TestMhTes:
    def test_mh_test(self):
        # self.assertEqual(expected, mh_test())
        raise SkipTest  # TODO: implement your test here


class TestBasicfit:
    def setUp(self):
        self.s1 = np.recarray(3, formats=['f8'], names=['p'])
        self.s2 = {'p': np.array([1, 2, 3])}

    def test_basicfit(self):
        assert_almost_equal(4.6666, Melding.basicfit(self.s1, self.s2), 1)


class TestMhTest:
    def test_mh_test(self):
        # assert_equal(expected, mh_test())
        raise SkipTest  # TODO: implement your test here


class TestClearNaN:
    def test_clear_na_n(self):
        obs = np.arange(9, dtype=float)
        obs.shape = (3, 3)
        obs[0, 2] = np.nan
        res = Melding.clearNaN(obs)
        print(res,  obs)
        assert_equal(res[0, 2],  np.mean([obs[0, 0], obs[0, 1]]) )
        

if __name__ == '__main__':
    unittest.main()
