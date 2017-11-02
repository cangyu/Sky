import unittest
from copy import deepcopy

import numpy as np
from src.iges import IGES_Model

from nurbs import BilinearSurf, Surf

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


def surf_insert_cmp(s, fn):
    """
    构造对比
    :param s: 原始曲面
    :type s: Surf
    :param fn: 文件名
    :return: None
    """

    ss = deepcopy(s)
    s.insert_knot(0.5, 1, 'U')
    s.insert_knot(0.5, 1, 'V')
    s.pan((0, 0, 10))

    model_file = IGES_Model(fn)
    model_file.add_entity(s.to_iges())
    model_file.add_entity(ss.to_iges())
    model_file.write()
    if auto_view:
        view(fn)

    print('\nOriginal:')
    print(ss)

    print('Inserted:')
    print(s)


class InsertTest(unittest.TestCase):
    @staticmethod
    def test():
        surf_insert_cmp(BilinearSurf(P1), 'insert_surf1.igs')
        surf_insert_cmp(BilinearSurf(P2), 'insert_surf2.igs')


if __name__ == '__main__':
    unittest.main()
