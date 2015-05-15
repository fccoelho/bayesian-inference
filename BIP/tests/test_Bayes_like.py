from __future__ import absolute_import
import unittest
from nose import SkipTest
from nose.tools import assert_equal
from numpy import arange,  sqrt
from numpy.testing import *
from numpy.random import lognormal
from scipy import stats
from BIP.Bayes.like import Normal, Lognormal, Poisson, Negbin, Binomial, Gamma, Weibull, Bernoulli, Beta, Categor

class TestCategor(unittest.TestCase):
    def test_categor(self):
        raise SkipTest # TODO: implement your test here

class TestNormal(unittest.TestCase):
    def test_normal(self):
        data = stats.norm.rvs(loc=0, scale=1, size=1000)
        ll = [Normal(data, i, 1) for i in arange(-15, 15, .001)]
        mle = arange(-15, 15, .001)[ll.index(max(ll))]
        #print "Normal MLE error: ", abs(0-mle)
        assert_almost_equal(mle, 0 ,1,"MLE, real: ",True)
        
class TestLognormal(unittest.TestCase):
    def test_lognormal(self):
        data = lognormal(5, 1, size=2000)
        ll = [Lognormal(data, i, 1) for i in arange(0.1, 15, .001)]
        mle = arange(0.1, 15, .001)[ll.index(max(ll))]
        #print "LogNormal MLE error: ", abs(5-mle)
        assert_almost_equal(mle, 5 ,1,"MLE, real: ",True)
        #assert abs(5-mle) <0.01

class TestPoisson(unittest.TestCase):
    def test_poisson(self):
        data = stats.poisson.rvs(5, size=2000)
        ll = [Poisson(data, i) for i in arange(0.1, 15, .001)]
        mle = arange(0.1, 15, .001)[ll.index(max(ll))]
        #print "Poisson MLE error: ", abs(5-mle)
        assert_almost_equal(mle, 5 ,1,"MLE, real: ",True)
        #assert abs(5-mle) <0.1

class TestNegbin(unittest.TestCase):
    def test_negbin(self):
        data = stats.nbinom .rvs(10, .5, size=1000)
        ll = [Negbin(data, 10, i) for i in arange(0.01, 1, .005)]
        mle = arange(0.01, 1, .005)[ll.index(max(ll))]
        #print "Negative Binomial MLE error: ", abs(.5-mle)
        assert_almost_equal(mle, .5 ,1,"MLE, real: ",True)
        #assert abs(.5-mle) <0.01

class TestBinomial(unittest.TestCase):
    def test_binomial(self):
        data = stats.binom.rvs(10, .5, size=1000)
        ll = [Binomial(data, 10, i) for i in arange(0.01, 1, .005)]
        mle = arange(0.01, 1, .005)[ll.index(max(ll))]
        #print "Binomial MLE error: ", abs(.5-mle)
        assert_almost_equal(mle, .5 ,1,"MLE, real: ",True)
        #assert abs(.5-mle) <0.01

class TestWeibull(unittest.TestCase):
    def test_weibull(self):
        data = stats.weibull_max(2, .5)
        #TODO:  find out which weibull function to use
        assert True 

class TestBernoulli(unittest.TestCase):
    def test_bernoulli(self):
        data = stats.bernoulli.rvs(.5, size=1000)
        ll = [Bernoulli(data, i) for i in arange(.01, 1, .001)]
        mle = arange(0.01, 1, .001)[ll.index(max(ll))]
        #print "Bernoulli MLE error: ", abs(.5-mle)
        assert_almost_equal(mle, .5 ,1,"MLE, real: ",True)
        #assert abs(.5-mle) <0.01

class TestGamma(unittest.TestCase):
    def test_gamma(self):
        data = stats.gamma.rvs(2,scale=2,  size=1000)
        ll = [Gamma(data, 2, i) for i in arange(.01, 5, .001)]
        mle = arange(0.01, 5, .001)[ll.index(max(ll))]
        #print "Gamma MLE error: ", abs(2-mle)
        assert_almost_equal(mle, 2 ,1,"MLE, real: ",True)
        #assert abs(2-mle) <0.01
class TestBeta(unittest.TestCase):
    def test_beta(self):
        data = stats.beta.rvs(2,2,  size=1000)
        ll = [Beta(data, 2, i) for i in arange(.01, 1, .001)]
        mle = arange(0.01, 1, .001)[ll.index(max(ll))]
        #print "Beta MLE error: ", abs(2-mle)
        assert_almost_equal(mle, 2 ,1,"MLE, real: ",True)
        #assert abs(2-mle) <0.01

class TestSimple(unittest.TestCase):
    def test_simple(self):
        raise SkipTest # TODO: implement your test here

if __name__ == '__main__':
    unittest.main()
