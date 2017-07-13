import unittest
import numpy as np
from src.nurbs.surface import BilinearSurf, ClampedNURBSSurf
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


def surf_swap_cmp(s, fn):
    """
    构造对比
    :param s: 原始曲面
    :type s: ClampedNURBSSurf
    :param fn: 文件名
    :return: None
    """

    ss = deepcopy(s)
    ss.elevate(1, 2)
    s.elevate(1, 2)
    s.swap()

    model_file = IGES_Model(fn)
    model_file.add_entity(s.to_iges())
    model_file.add_entity(ss.to_iges())
    model_file.write()
    if auto_view:
        view(fn)

    print('\nBefore:')
    print(ss)

    print('After:')
    print(s)

    I, J = 50, 30
    u_dist, v_dist = np.meshgrid(np.linspace(0, 1, I), np.linspace(0, 1, J), indexing='ij')
    ps = np.zeros((I, J, 3))
    pss = np.zeros((I, J, 3))

    for i in range(I):
        for j in range(J):
            ps[i][j] = s(v_dist[i][j], u_dist[i][j])
            pss[i][j] = ss(u_dist[i][j], v_dist[i][j])

    print('\n')
    print(np.allclose(ps, pss))


class SwapTest(unittest.TestCase):
    @staticmethod
    def test():
        surf_swap_cmp(BilinearSurf(P1), 'swap_surf1.igs')
        surf_swap_cmp(BilinearSurf(P2), 'swap_surf2.igs')


if __name__ == '__main__':
    unittest.main()
