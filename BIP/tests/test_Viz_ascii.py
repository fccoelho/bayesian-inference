from __future__ import absolute_import
from nose import SkipTest
from nose.tools import assert_equal
from BIP.Viz.ascii import Histogram

from numpy.random import normal,  seed
seed(1)
data = normal(size=1000)
bins = 10
rnge = (-1, 1)
height = 10
character = '|'

class TestHistogram:
    def test___init__(self):
        histogram = Histogram(data, bins, rnge)
        raise SkipTest # TODO: implement your test here

    def test_horizontal(self):
        histogram = Histogram(data, bins, rnge)
        assert_equal("1.00\n", histogram.horizontal(height, character)[-5:])
        

    def test_vertical(self):
        histogram = Histogram(data, bins, rnge)
        assert_equal("0.80", histogram.vertical(height, character).splitlines()[-1][:4])
