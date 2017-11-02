import unittest
from copy import deepcopy

import numpy as np

from nurbs import BilinearSurf, Surf
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


def surf_pan_cmp(s, dv, fn):
    """
    构造对比
    :param s: 原始曲面
    :type s: Surf
    :param dv: 偏移向量
    :param fn: 文件名
    :return: None
    """

    ss = deepcopy(s)
    ss.pan(dv)

    model_file = IGES_Model(fn)
    model_file.add_entity(s.to_iges())
    model_file.add_entity(ss.to_iges())
    model_file.write()
    if auto_view:
        view(fn)


class PanTest(unittest.TestCase):
    @staticmethod
    def test():
        surf_pan_cmp(BilinearSurf(P1), (0, 0, 10), 'pan_surf1.igs')
        surf_pan_cmp(BilinearSurf(P2), (0, 0, 10), 'pan_surf2.igs')


if __name__ == '__main__':
    unittest.main()
