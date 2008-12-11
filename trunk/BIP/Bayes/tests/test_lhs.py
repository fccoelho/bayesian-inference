import unittest
from BIP.Bayes import lhs
from scipy import stats

class TestLhs(unittest.TestCase):
    def test_lhs(self):
        s = lhs.lhs(stats.norm,(0,1),size=(10,100))
        assert s.shape == (10,100)

if __name__ == '__main__':
    unittest.main()
