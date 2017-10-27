import unittest
from copy import deepcopy

import numpy as np

from src.iges import IGES_Model
from src.nurbs.surface import BilinearSurf, ClampedNURBSSurf

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


def surf_refine_cmp(s, fn):
    """
    构造对比
    :param s: 原始曲面
    :type s: ClampedNURBSSurf
    :param fn: 文件名
    :return: None
    """

    ss = deepcopy(s)
    s.refine('U', [0.2, 0.3, 0.6, 0.7])
    s.refine('V', [0.1, 0.4, 0.8, 0.9])

    model_file = IGES_Model(fn)
    model_file.add_entity(s.to_iges())
    model_file.add_entity(ss.to_iges())
    model_file.write()
    if auto_view:
        view(fn)

    print('\nOriginal:')
    print(ss)

    print('Refined:')
    print(s)

    I, J = 50, 30
    u_dist, v_dist = np.meshgrid(np.linspace(0, 1, I), np.linspace(0, 1, J), indexing='ij')
    ps = np.zeros((I, J, 3))
    pss = np.zeros((I, J, 3))

    for i in range(I):
        for j in range(J):
            ps[i][j] = s(u_dist[i][j], v_dist[i][j])
            pss[i][j] = ss(u_dist[i][j], v_dist[i][j])

    print('\n')
    print(np.allclose(ps, pss))


class RefineTest(unittest.TestCase):
    @staticmethod
    def test():
        surf_refine_cmp(BilinearSurf(P1), 'refine_surf1.igs')
        surf_refine_cmp(BilinearSurf(P2), 'refine_surf2.igs')


if __name__ == '__main__':
    unittest.main()
