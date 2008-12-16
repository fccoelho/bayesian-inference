import unittest
from BIP.Bayes import lhs
from scipy import stats

class TestLhs(unittest.TestCase):
    def test_lhs(self):
        s = lhs.lhs(stats.norm,(0,1),size=(10,100))
        assert s.shape == (10,100)

class TestLhsFromSample(unittest.TestCase):
    def test_lhs_from_sample(self):
        assert False # TODO: implement your test here

class TestLhsFromDensity(unittest.TestCase):
    def test_lhs_from_density(self):
        assert False # TODO: implement your test here

if __name__ == '__main__':
    unittest.main()
