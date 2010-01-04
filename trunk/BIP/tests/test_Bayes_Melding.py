import unittest


class TestFitModel(unittest.TestCase):
    def test___init__(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        assert False # TODO: implement your test here

    def test_do_inference(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.do_inference(prior, data, predlen, method))
        assert False # TODO: implement your test here

    def test_init_priors(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.init_priors(prior))
        assert False # TODO: implement your test here

    def test_monitor(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.monitor())
        assert False # TODO: implement your test here

    def test_plot(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.plot())
        assert False # TODO: implement your test here

    def test_run(self):
        # fit_model = FitModel(K, L, model, ntheta, nphi, inits, tf, phinames, thetanames, verbose)
        # self.assertEqual(expected, fit_model.run(wl, nw, data, method, monitor))
        assert False # TODO: implement your test here

class TestMeld(unittest.TestCase):
    def test___init__(self):
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        assert False # TODO: implement your test here

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
        # meld = Meld(K, L, model, ntheta, nphi, alpha, verbose, viz)
        # self.assertEqual(expected, meld.sir(data, t, variance, nopool, savetemp))
        assert False # TODO: implement your test here

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