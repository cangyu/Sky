import unittest
import numpy as np
import time
from src.nurbs.surface import BilinearSurf
from src.iges.iges_core import IGES_Model

L = 10
P1 = np.array([[[0, 0, 0], [0, L, L]],
               [[L, 0, L], [L, L, 0]]], float)

P2 = np.array([[[0, 0, L], [0, L, 0]],
               [[L, 0, L], [L, L, 0]]], float)


def build_bilinear_surf(P):
    bsf = BilinearSurf(P)
    suffix = time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())
    model_file = IGES_Model('bilinear_{}.igs'.format(suffix))
    model_file.add_entity(bsf.to_iges())
    model_file.write()


class nurbs_surf_basic_test(unittest.TestCase):
    @staticmethod
    def test_bilinear():
        # build_bilinear_surf(P1)
        build_bilinear_surf(P2)
