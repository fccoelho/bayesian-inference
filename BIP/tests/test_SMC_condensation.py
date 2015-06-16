from __future__ import absolute_import
import unittest
from nose import SkipTest


class TestGlobalData(unittest.TestCase):
    def test___init__(self):
        # global_data = GlobalData(prior, process, scene, observ, nsam, nit)
        raise SkipTest # TODO: implement your test here

class TestIterationData(unittest.TestCase):
    def test___init__(self):
        # iteration_data = IterationData()
        raise SkipTest # TODO: implement your test here

class TestCondensation(unittest.TestCase):
    def test___init__(self):
        # condensation = Condensation(model)
        raise SkipTest # TODO: implement your test here

    def test_calculateBaseWeights(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.calculateBaseWeights())
        raise SkipTest # TODO: implement your test here

    def test_pickBaseSample(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.pickBaseSample())
        raise SkipTest # TODO: implement your test here

    def test_predictNewBases(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.predictNewBases())
        raise SkipTest # TODO: implement your test here

    def test_runFilter(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.runFilter())
        raise SkipTest # TODO: implement your test here

    def test_updateAfterIterating(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.updateAfterIterating(iteration))
        raise SkipTest # TODO: implement your test here

class TestEvaluateGaussian(unittest.TestCase):
    def test_evaluate_gaussian(self):
        # self.assertEqual(expected, evaluate_gaussian(val, sigma))
        raise SkipTest # TODO: implement your test here

if __name__ == '__main__':
    unittest.main()
