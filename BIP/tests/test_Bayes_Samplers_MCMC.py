from __future__ import absolute_import
from nose import SkipTest
from nose.tools import assert_equal
from numpy.testing import *
from BIP.Bayes.Samplers.MCMC import *
import unittest
 

class TestMetropolis:
    def test___init__(self):
        # metropolis = Metropolis(proposal_dist, likfun)
        raise SkipTest # TODO: implement your test here

    def test_add_salt(self):
        # metropolis = Metropolis(proposal_dist, likfun)
        # assert_equal(expected, metropolis.add_salt(dataset, band))
        raise SkipTest # TODO: implement your test here

    def test_propose(self):
        # metropolis = Metropolis(proposal_dist, likfun)
        # assert_equal(expected, metropolis.propose())
        raise SkipTest # TODO: implement your test here

    def test_step(self):
        # metropolis = Metropolis(proposal_dist, likfun)
        # assert_equal(expected, metropolis.step(n))
        raise SkipTest # TODO: implement your test here

    def test_tune(self):
        # metropolis = Metropolis(proposal_dist, likfun)
        # assert_equal(expected, metropolis.tune(ar))
        raise SkipTest # TODO: implement your test here

class test__Sampler:
    def TestDIC(self):
        # __sampler = _Sampler(parpriors, parnames)
        # assert_equal(expected, __sampler.DIC())
        raise SkipTest # TODO: implement your test here

    def test___init__(self):
        # __sampler = _Sampler(parpriors, parnames)
        raise SkipTest # TODO: implement your test here

    def test_dimensions(self):
        # __sampler = _Sampler(parpriors, parnames)
        # assert_equal(expected, __sampler.dimensions())
        raise SkipTest # TODO: implement your test here

    def test_gr_R(self):
        # __sampler = _Sampler(parpriors, parnames)
        # assert_equal(expected, __sampler.gr_R(end, start))
        raise SkipTest # TODO: implement your test here

    def test_gr_convergence(self):
        # __sampler = _Sampler(parpriors, parnames)
        # assert_equal(expected, __sampler.gr_convergence(relevantHistoryEnd, relevantHistoryStart))
        raise SkipTest # TODO: implement your test here

    def test_po(self):
        # __sampler = _Sampler(parpriors, parnames)
        # assert_equal(expected, __sampler.po())
        raise SkipTest # TODO: implement your test here

    def test_setup_xmlrpc_plotserver(self):
        # __sampler = _Sampler(parpriors, parnames)
        # assert_equal(expected, __sampler.setup_xmlrpc_plotserver())
        raise SkipTest # TODO: implement your test here

    def test_term_pool(self):
        # __sampler = _Sampler(parpriors, parnames)
        # assert_equal(expected, __sampler.term_pool())
        raise SkipTest # TODO: implement your test here

class TestModelAsRa:
    def test_model_as_ra(self):
        # assert_equal(expected, model_as_ra(theta, model, phinames))
        raise SkipTest # TODO: implement your test here

class TestDream:
    def test___init__(self):
        # dream = Dream(meldobj, samples, sampmax, data, t, parpriors, parnames, parlimits, likfun, likvariance, burnin, thin, convergenceCriteria, nCR, DEpairs, adaptationRate, eps, mConvergence, mAccept, **kwargs)
        raise SkipTest # TODO: implement your test here

    def test_delayed_rejection(self):
        # dream = Dream(meldobj, samples, sampmax, data, t, parpriors, parnames, parlimits, likfun, likvariance, burnin, thin, convergenceCriteria, nCR, DEpairs, adaptationRate, eps, mConvergence, mAccept, **kwargs)
        # assert_equal(expected, dream.delayed_rejection(xi, zi, pxi, zprob))
        raise SkipTest # TODO: implement your test here

    def test_step(self):
        # dream = Dream(meldobj, samples, sampmax, data, t, parpriors, parnames, parlimits, likfun, likvariance, burnin, thin, convergenceCriteria, nCR, DEpairs, adaptationRate, eps, mConvergence, mAccept, **kwargs)
        # assert_equal(expected, dream.step())
        raise SkipTest # TODO: implement your test here

    def test_update_CR_dist(self):
        # dream = Dream(meldobj, samples, sampmax, data, t, parpriors, parnames, parlimits, likfun, likvariance, burnin, thin, convergenceCriteria, nCR, DEpairs, adaptationRate, eps, mConvergence, mAccept, **kwargs)
        # assert_equal(expected, dream.update_CR_dist())
        raise SkipTest # TODO: implement your test here

