import unittest
from numpy import arange
from numpy.random import lognormal
from scipy import stats
from Bayes.like import Normal, Lognormal, Poisson, Negbin, Binomial, Gamma, Weibull, Bernoulli, Beta, Categor

class TestCategor(unittest.TestCase):
    def test_categor(self):
        assert False # TODO: implement your test here

class TestNormal(unittest.TestCase):
    def test_normal(self):
        data = stats.norm.rvs(loc=0, scale=1, size=1000)
        ll = [Normal(data, i, 1) for i in arange(-15, 15, .005)]
        mle = arange(-15, 15, .005)[ll.index(max(ll))]
        print "Normal MLE error: ", abs(0-mle)
        assert abs(0-mle) <0.05
class TestLognormal(unittest.TestCase):
    def test_lognormal(self):
        data = lognormal(5, 1, size=1000)
        ll = [Lognormal(data, i, 1) for i in arange(0.1, 15, .005)]
        mle = arange(0.1, 15, .005)[ll.index(max(ll))]
        print "LogNormal MLE error: ", abs(5-mle)
        assert abs(5-mle) <0.01

class TestPoisson(unittest.TestCase):
    def test_poisson(self):
        data = stats.poisson.rvs(5, size=1000)
        ll = [Poisson(data, i) for i in arange(0.1, 15, .005)]
        mle = arange(0.1, 15, .005)[ll.index(max(ll))]
        print "Poisson MLE error: ", abs(5-mle)
        assert abs(5-mle) <0.05

class TestNegbin(unittest.TestCase):
    def test_negbin(self):
        data = stats.nbinom .rvs(10, .5, size=1000)
        ll = [Negbin(data, 10, i) for i in arange(0.01, 1, .005)]
        mle = arange(0.01, 1, .005)[ll.index(max(ll))]
        print "Negative Binomial MLE error: ", abs(.5-mle)
        assert abs(.5-mle) <0.01

class TestBinomial(unittest.TestCase):
    def test_binomial(self):
        data = stats.binom.rvs(10, .5, size=1000)
        ll = [Binomial(data, 10, i) for i in arange(0.01, 1, .005)]
        mle = arange(0.01, 1, .005)[ll.index(max(ll))]
        print "Binomial MLE error: ", abs(.5-mle)
        assert abs(.5-mle) <0.01

class TestWeibull(unittest.TestCase):
    def test_weibull(self):
        data = stats.weibull_max(2, .5)
        #TODO:  find out which weibull function to use
        assert True 

class TestBernoulli(unittest.TestCase):
    def test_bernoulli(self):
        assert False # TODO: implement your test here

class TestGamma(unittest.TestCase):
    def test_gamma(self):
        assert False # TODO: implement your test here

class TestBeta(unittest.TestCase):
    def test_beta(self):
        assert False # TODO: implement your test here

class TestSimple(unittest.TestCase):
    def test_simple(self):
        assert False # TODO: implement your test here

if __name__ == '__main__':
    unittest.main()
