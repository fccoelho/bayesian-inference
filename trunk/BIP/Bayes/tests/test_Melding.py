import unittest
from BIP.Bayes import Melding

class TestMeld(unittest.TestCase):
    def test_abcRun(self):
        assert False # TODO: implement your test here

    def test_addData(self):
        assert False # TODO: implement your test here

    def test_basicfit(self):
        assert False # TODO: implement your test here

    def test_filtM(self):
        assert False # TODO: implement your test here

    def test_getPosteriors(self):
        assert False # TODO: implement your test here

    def test_logPooling(self):
        assert False # TODO: implement your test here

    def test_object_initialization(self):
        Me = Melding.Meld(K=1000,L=1000,model=lambda x:x, ntheta=4,nphi=5)
        assert isinstance(Me, Melding.Meld)

    def test_run(self):
        assert False # TODO: implement your test here

    def test_setPhi(self):
        assert False # TODO: implement your test here

    def test_setPhiFromData(self):
        assert False # TODO: implement your test here

    def test_setTheta(self):
        assert False # TODO: implement your test here

    def test_setThetaFromData(self):
        assert False # TODO: implement your test here

    def test_sir(self):
        assert False # TODO: implement your test here

class TestModel(unittest.TestCase):
    def test_model(self):
        assert False # TODO: implement your test here

class TestRun(unittest.TestCase):
    def test_run(self):
        assert False # TODO: implement your test here

class TestKDE(unittest.TestCase):
    def test_k_d_e(self):
        assert False # TODO: implement your test here

class TestLikeli(unittest.TestCase):
    def test_likeli(self):
        assert False # TODO: implement your test here

class TestFilt(unittest.TestCase):
    def test_filt(self):
        assert False # TODO: implement your test here

class TestFiltM(unittest.TestCase):
    def test_filt_m(self):
        assert False # TODO: implement your test here

class TestSIR(unittest.TestCase):
    def test_s_i_r(self):
        assert False # TODO: implement your test here

class TestMain(unittest.TestCase):
    def test_main(self):
        assert False # TODO: implement your test here

class TestMain2(unittest.TestCase):
    def test_main2(self):
        assert False # TODO: implement your test here

if __name__ == '__main__':
    unittest.main()
