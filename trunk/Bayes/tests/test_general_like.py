from nose import SkipTest
from scipy import stats
from numpy import arange
#from Bayes.general.like import *

class TestCategor:
    def test_categor(self):
        raise SkipTest # TODO: implement your test here

class TestNormal:
    def test_normal(self):
        data = norm(loc=0, scale=1, size=100)
        mle = max([Normal(data, i, 1) for i in arange(-2, 2, .1)])
        print mle


class TestLognormal:
    def test_lognormal(self):
        raise SkipTest # TODO: implement your test here

class TestPoisson:
    def test_poisson(self):
        raise SkipTest # TODO: implement your test here

class TestNegbin:
    def test_negbin(self):
        raise SkipTest # TODO: implement your test here

class TestBinomial:
    def test_binomial(self):
        raise SkipTest # TODO: implement your test here

class TestWeibull:
    def test_weibull(self):
        raise SkipTest # TODO: implement your test here

class TestBernoulli:
    def test_bernoulli(self):
        raise SkipTest # TODO: implement your test here

class TestGamma:
    def test_gamma(self):
        raise SkipTest # TODO: implement your test here

class TestBeta:
    def test_beta(self):
        raise SkipTest # TODO: implement your test here

class TestSimple:
    def test_simple(self):
        raise SkipTest # TODO: implement your test here

