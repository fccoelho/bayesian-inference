import unittest
from numpy import arange
from scipy import stats

class TestCategor(unittest.TestCase):
    def test_categor(self):
        assert False # TODO: implement your test here

class TestNormal(unittest.TestCase):
    def test_normal(self):
        data = stats.norm.rvs(loc=0, scale=1, size=100)
        mle = [Normal(data, i, 1) for i in arange(-2, 2, .1)]
        mle = arange(-2, 2, .1)[mle.index(max(mle))]
        assert 0-mle <0.000001
class TestLognormal(unittest.TestCase):
    def test_lognormal(self):
        assert False # TODO: implement your test here

class TestPoisson(unittest.TestCase):
    def test_poisson(self):
        assert False # TODO: implement your test here

class TestNegbin(unittest.TestCase):
    def test_negbin(self):
        assert False # TODO: implement your test here

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
