import unittest
import numpy as np
from src.nurbs.nurbs_surface import NURBS_Surface
from src.iges.iges_core import IGES_Model


class nurbs_surf_basic_test(unittest.TestCase):
    def test_bilinear(self):
        L = 10
        U = np.array([0, 0, 1, 1])
        V = np.array([0, 0, 1, 1])
        Pw = np.array([[[L, 0, L, 1], [0, 0, 0, 1]], [[L, L, 0, 1], [0, L, L, 1]]])
        bsf = NURBS_Surface(U, V, Pw)

        model_file = IGES_Model('bilinear.igs')
        model_file.add_entity(bsf.to_iges(0, 0, 0, 0))
        model_file.write()
