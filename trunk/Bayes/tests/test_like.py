import unittest
from numpy import arange
from numpy.random import lognormal
from scipy import stats
from Bayes.like import Normal, Lognormal, Poisson

class TestCategor(unittest.TestCase):
    def test_categor(self):
        assert False # TODO: implement your test here

class TestNormal(unittest.TestCase):
    def test_normal(self):
        data = stats.norm.rvs(loc=0, scale=1, size=1000)
        mle = [Normal(data, i, 1) for i in arange(-2, 2, .1)]
        mle = arange(-2, 2, .1)[mle.index(max(mle))]
        print "Normal MLE error: ", (0-mle)
        assert abs(0-mle) <0.000001
class TestLognormal(unittest.TestCase):
    def test_lognormal(self):
        data = lognormal.rvs(loc=5, scale=1, size=1000)
        mle = [Poisson(data, i) for i in arange(0, 10, .1)]
        mle = arange(0, 10, .1)[mle.index(max(mle))]
        print "LogNormal MLE error: ", (5-mle)
        assert abs(5-mle) <0.000001

class TestPoisson(unittest.TestCase):
    def test_poisson(self):
        data = stats.poisson.rvs(5, size=1000)
        ll = [Poisson(data, i) for i in arange(0, 10, .1)]
        mle = arange(0, 10, .1)[ll.index(max(ll))]
        print "Poisson MLE error: ", (5-mle)
        assert abs(5-mle) <0.000001

class TestNegbin(unittest.TestCase):
    def test_negbin(self):
        #TODO Unfinished! 
        data = stats.nbinom .rvs(n=10, pr=.5, size=1000)
        mle = [Negbin(data, 10, i) for i in arange(0, 1, .01)]
        mle = arange(0, 1, .01)[mle.index(max(mle))]
        print "Negative Binomial MLE error: ", (.5-mle)
        assert abs(.5-mle) <0.000001

class TestBinomial(unittest.TestCase):
    def test_binomial(self):
        assert False # TODO: implement your test here

class TestWeibull(unittest.TestCase):
    def test_weibull(self):
        assert False # TODO: implement your test here

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
