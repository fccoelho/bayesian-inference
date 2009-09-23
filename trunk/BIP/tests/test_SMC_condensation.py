import unittest

class TestGlobalData(unittest.TestCase):
    def test___init__(self):
        # global_data = GlobalData(prior, process, scene, observ, nsam, nit)
        assert False # TODO: implement your test here

class TestIterationData(unittest.TestCase):
    def test___init__(self):
        # iteration_data = IterationData()
        assert False # TODO: implement your test here

class TestCondensation(unittest.TestCase):
    def test___init__(self):
        # condensation = Condensation(model)
        assert False # TODO: implement your test here

    def test_calculateBaseWeights(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.calculateBaseWeights())
        assert False # TODO: implement your test here

    def test_pickBaseSample(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.pickBaseSample())
        assert False # TODO: implement your test here

    def test_predictNewBases(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.predictNewBases())
        assert False # TODO: implement your test here

    def test_runFilter(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.runFilter())
        assert False # TODO: implement your test here

    def test_updateAfterIterating(self):
        # condensation = Condensation(model)
        # self.assertEqual(expected, condensation.updateAfterIterating(iteration))
        assert False # TODO: implement your test here

class TestEvaluateGaussian(unittest.TestCase):
    def test_evaluate_gaussian(self):
        # self.assertEqual(expected, evaluate_gaussian(val, sigma))
        assert False # TODO: implement your test here

if __name__ == '__main__':
    unittest.main()
