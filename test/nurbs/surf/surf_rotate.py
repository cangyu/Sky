import unittest
import numpy as np
from src.nurbs.surface import BilinearSurf
from src.iges.iges_core import IGES_Model
from copy import deepcopy

try:
    from src.misc.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True

L = 10
P1 = np.array([[[0, 0, 0], [0, L, L]],
               [[L, 0, L], [L, L, 0]]], float)

P2 = np.array([[[0, 0, L], [0, L, 0]],
               [[L, 0, L], [L, L, 0]]], float)


def surf_rotate_cmp(s, ref, ax, ang, fn):
    ss = deepcopy(s)
    ss.rotate(ref, ax, ang)

    model_file = IGES_Model(fn)
    model_file.add_entity(s.to_iges())
    model_file.add_entity(ss.to_iges())
    model_file.write()
    if auto_view:
        view(fn)


class RotateTest(unittest.TestCase):
    @staticmethod
    def test():
        surf_rotate_cmp(BilinearSurf(P1), [L, L, L], [0, 0, 5], 45, 'rotate_surf1.igs')
        surf_rotate_cmp(BilinearSurf(P2), [L, L, L], [0, 0, 5], 45, 'rotate_surf2.igs')


if __name__ == '__main__':
    unittest.main()
