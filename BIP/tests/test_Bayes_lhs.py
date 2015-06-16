from __future__ import absolute_import
import unittest
from nose import SkipTest
from nose.tools import assert_equal
from numpy.testing import *
from BIP.Bayes import like, lhs
import scipy.stats as st
import numpy


class TestLhsFromSample(unittest.TestCase):
    def test_lhs_from_sample(self):
        # self.assertEqual(expected, lhsFromSample(sample, siz))
        raise SkipTest  # TODO: implement your test here


class TestLhsFromDensity(unittest.TestCase):
    def test_lhs_from_density(self):
        # self.assertEqual(expected, lhsFromDensity(kde, siz))
        raise SkipTest  # TODO: implement your test here


class TestLhs(unittest.TestCase):
    def test_lhs(self):
        # self.assertEqual(expected, lhs(dist, parms, siz))
        raise SkipTest  # TODO: implement your test here

    def test_correlation_structure(self):
        cm = numpy.array([[1, .8], [.8, 1]])
        c = lhs.lhs([st.norm, st.beta], [(50, 1), (10, 2)], 2000, False, cm)
        r, p = st.spearmanr(c[0], c[1])
        assert_almost_equal(r, 0.8, 1)
        assert p < 0.05

    def test_adherence_to_target_dist(self):
        raise SkipTest


if __name__ == '__main__':
    unittest.main()
