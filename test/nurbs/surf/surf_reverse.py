import unittest
from copy import deepcopy

import numpy as np
from src.iges import IGES_Model

from nurbs import BilinearSurf

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


def surf_reverse_cmp(s, d, fn):
    ss = deepcopy(s)
    ss.reverse(d)

    model_file = IGES_Model(fn)
    model_file.add_entity(s.to_iges())
    model_file.add_entity(ss.to_iges())
    model_file.write()
    if auto_view:
        view(fn)

    print('Origin:')
    print(s)

    print('Reversed:')
    print(ss)


class ReverseTest(unittest.TestCase):
    @staticmethod
    def test():
        surf_reverse_cmp(BilinearSurf(P1), 'U', 'reverse_surf1_U.igs')
        surf_reverse_cmp(BilinearSurf(P1), 'V', 'reverse_surf1_V.igs')
        surf_reverse_cmp(BilinearSurf(P1), 'UV', 'reverse_surf1_UV.igs')

        surf_reverse_cmp(BilinearSurf(P2), 'U', 'reverse_surf2_U.igs')
        surf_reverse_cmp(BilinearSurf(P2), 'V', 'reverse_surf2_V.igs')
        surf_reverse_cmp(BilinearSurf(P2), 'UV', 'reverse_surf2_UV.igs')


if __name__ == '__main__':
    unittest.main()
