import unittest

import numpy as np

from nurbs import BilinearSurf
from src.iges import IGES_Model

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

L = 10
P1 = np.array([[[0, 0, 0], [0, L, L]],
               [[L, 0, L], [L, L, 0]]], float)

P2 = np.array([[[0, 0, L], [0, L, 0]],
               [[L, 0, L], [L, L, 0]]], float)


def build_bilinear_surf(P, fn):
    bsf = BilinearSurf(P)
    model_file = IGES_Model(fn)
    model_file.add_entity(bsf.to_iges(0, 0, 0, 0))
    model_file.write()
    view(fn)


class BasicSurfTest(unittest.TestCase):
    @staticmethod
    def test_bilinear():
        build_bilinear_surf(P1, 'surf1.igs')
        build_bilinear_surf(P2, 'surf2.igs')


if __name__ == '__main__':
    unittest.main()
