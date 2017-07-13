import unittest
import numpy as np
from src.nurbs.surface import BilinearSurf, ClampedNURBSSurf
from src.iges.iges_core import IGES_Model
from copy import deepcopy
from src.nurbs.utility import equal
from numpy.linalg import norm

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


def surf_extract_cmp(s, fn):
    """
    构造对比
    :param s: 原始曲面
    :type s: ClampedNURBSSurf
    :param fn: 文件名
    :return: None
    """

    mid = s.extract('U', 0.5)

    model_file = IGES_Model(fn)
    model_file.add_entity(s.to_iges())
    model_file.add_entity(mid.to_iges())
    model_file.write()
    if auto_view:
        view(fn)

    N = 1000
    vl = np.linspace(0, 1, N)
    vvals = np.zeros((N, 3))
    vvall = np.zeros((N, 3))

    for i in range(N):
        vvals[i] = s(0.5, vl[i])
        vvall[i] = mid(vl[i])

    print(np.allclose(vvals, vvall))


class ExtractTest(unittest.TestCase):
    @staticmethod
    def test():
        surf_extract_cmp(BilinearSurf(P1), 'extract_surf1.igs')
        surf_extract_cmp(BilinearSurf(P2), 'extract_surf2.igs')


if __name__ == '__main__':
    unittest.main()
