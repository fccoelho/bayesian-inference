from __future__ import absolute_import
import unittest
from BIP.Viz.realtime import RTplot
from numpy import *

data = random.normal(size=(10,100))
title = 'test'
names = []
style = 'lines'

class TestRTplot(unittest.TestCase):
    def test___init__(self):
        r_tplot = RTplot(persist=0)

    def test_clearFig(self):
        r_tplot = RTplot(persist=0)
        self.assertEqual(None, r_tplot.clearFig())
        

    def test_plothist(self):
        r_tplot = RTplot(persist=0)
        self.assertEqual(None, r_tplot.plothist(data, title, names))

    def test_plotlines(self):
        r_tplot = RTplot(persist=0)
        self.assertEqual(None, r_tplot.plotlines(data, names, title, style))
     

    def test_scatter(self):
        r_tplot = RTplot(persist=0)
        self.assertEqual(None, r_tplot.scatter([1,2,3], [1,2,3], title='test'))
     

if __name__ == '__main__':
    unittest.main()
