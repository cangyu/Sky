import unittest

import numpy as np

from nurbs import GlobalInterpolatedSurf
from src.iges import IGES_Model


def build_wing(foil, N: int, L: float, p: int, q: int):
    """
    构建直机翼
    :param foil: 翼型名称
    :param N: 段数量
    :param L: 每段长度
    :param p: U方向次数
    :param q: v方向次数
    :return: None
    """

    '''Read airfoil points'''
    fn = '../../airfoil/' + foil + '.dat'
    pnt_list = []
    fin = open(fn)
    for pnt in fin:
        x, y, z = pnt.split()
        pnt_list.append([float(x), float(y), float(z)])
    fin.close()
    pts = np.copy(pnt_list)

    '''Assemble all points'''
    N += 1
    m, dim = pts.shape
    all_pts = np.zeros((N, m, dim))
    for i in range(0, N):
        all_pts[i] = np.copy(pts)
        for j in range(0, m):
            all_pts[i][j][-1] = L * i

    '''Build IGS File'''
    wsf = GlobalInterpolatedSurf(all_pts, p, q)
    model_file = IGES_Model("{}_{}_{}_{}.igs".format(foil, p, q, (N - 1) * L))
    model_file.add_entity(wsf.to_iges(0, 0, 0, 0))
    model_file.write()


class SurfGlobalInterpolationTest(unittest.TestCase):
    @staticmethod
    def test_wing():
        build_wing('M6', 10, 0.7, 3, 5)
        build_wing('M6', 10, 0.7, 5, 5)
        build_wing('NACA0012', 10, 0.7, 3, 5)
        build_wing('NACA0012', 10, 0.7, 5, 5)
        build_wing('RAE2822', 10, 0.7, 3, 5)
        build_wing('RAE2822', 10, 0.7, 5, 5)
