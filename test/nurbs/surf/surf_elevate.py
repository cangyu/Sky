import unittest
import numpy as np
from src.nurbs.surface import BilinearSurf, ClampedNURBSSurf
from src.iges.iges_core import IGES_Model
from copy import deepcopy

try:
    from src.com.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True

L = 10
P1 = np.array([[[0, 0, 0], [0, L, L]],
               [[L, 0, L], [L, L, 0]]], float)

P2 = np.array([[[0, 0, L], [0, L, 0]],
               [[L, 0, L], [L, L, 0]]], float)


def surf_elevate_cmp(s, fn):
    """
    构造对比
    :param s: 原始曲面
    :type s: ClampedNURBSSurf
    :param fn: 文件名
    :return: None
    """

    ss = deepcopy(s)
    s.elevate(1, 2)

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
            ps[i][j] = s(u_dist[i][j], v_dist[i][j])
            pss[i][j] = ss(u_dist[i][j], v_dist[i][j])

    print('\n')
    print(np.allclose(ps, pss))


class ElevateTest(unittest.TestCase):
    @staticmethod
    def test():
        surf_elevate_cmp(BilinearSurf(P1), 'elevate_surf1.igs')
        surf_elevate_cmp(BilinearSurf(P2), 'elevate_surf2.igs')


if __name__ == '__main__':
    unittest.main()
