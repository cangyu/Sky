import unittest
import numpy as np
from src.nurbs.surface import BilinearSurf
from src.iges.iges_core import IGES_Model
from copy import deepcopy

try:
    from src.com.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = False

L = 10
P1 = np.array([[[0, 0, 0], [0, L, L]],
               [[L, 0, L], [L, L, 0]]], float)

P2 = np.array([[[0, 0, L], [0, L, 0]],
               [[L, 0, L], [L, L, 0]]], float)


def build_bilinear_surf(P, fn):
    bsf0 = BilinearSurf(P)
    bsf1 = deepcopy(bsf0)
    bsf1.rotate([0, 0, 0], [5, 5, 5], 45)
    model_file = IGES_Model(fn)
    model_file.add_entity(bsf0.to_iges())
    model_file.add_entity(bsf1.to_iges())
    model_file.write()
    if auto_view:
        view(fn)


class BasicSurfTest(unittest.TestCase):
    @staticmethod
    def test_bilinear():
        build_bilinear_surf(P1, 'surf1.igs')
        build_bilinear_surf(P2, 'surf2.igs')


if __name__ == '__main__':
    unittest.main()
