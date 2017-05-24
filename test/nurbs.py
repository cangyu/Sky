import unittest
import numpy as np
from src.nurbs.basis import Basis


class nurbs_test(unittest.TestCase):
    def test_basis(self):
        U = np.array([0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5])
        N = Basis(U, 0)
        self.assertEquals(N(0, 0, 1.5), 0)
        self.assertEquals(N(1, 0, 1.5), 0)
        self.assertEquals(N(3, 0, 1.5), 1)
        self.assertEquals(N(7, 0, 4), 1)
        self.assertEquals(N(7, 0, 5), 0)
        self.assertEquals(N(1, 1, 0.1), 0.9)
        self.assertEquals(N(7, 2, 4.5), 0.25)
        self.assertEquals(N(7, 2, 5), 0)
        self.assertEquals(N(4, 2, 2.5), 0.125)
        self.assertEquals(N(4, 2, 2.5, 1), 0.5)
        self.assertEquals(N(4, 2, 2.5, 2), 1)
        self.assertEquals(N(3, 2, 2.5), 0.75)

    def test_crv(self):
        pass
