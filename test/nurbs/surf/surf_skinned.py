import unittest

import numpy as np
from src.iges import IGES_Model

from nurbs import GlobalInterpolatedCrv
from nurbs import Skinned
from src.aircraft.wing import Airfoil

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True


def write_wing(foil, N, L, p, q, with_foil=True):
    """
    构造机翼
    :param foil:翼型名称
    :param N: 分段数
    :param L: 每段长度
    :param p: 翼型曲线次数
    :param q: 展向插值次数
    :return: None
    """

    N += 1
    pts = Airfoil.read_pts(foil)
    m, dim = pts.shape
    all_pts = np.zeros((N, m, dim))
    for i in range(N):
        all_pts[i] = np.copy(pts)
        for j in range(0, m):
            all_pts[i][j][-1] = L * i

    crv_list = []
    for i in range(N):
        crv_list.append(GlobalInterpolatedCrv(all_pts[i], p))

    wsf = Skinned(crv_list, p, q)
    fn = "{}_{}_{}_{}_Skinning.igs".format(foil, p, q, (N - 1) * L)
    model_file = IGES_Model(fn)
    model_file.add_entity(wsf.to_iges())

    if with_foil:
        for crv in crv_list:
            model_file.add_entity(crv.to_iges(1, 0, [0, 0, 1]))

    model_file.write()

    if auto_view:
        view(fn)


class SkinnedTest(unittest.TestCase):
    @staticmethod
    def test_wing():
        write_wing('M6', 10, 0.7, 3, 5)
        write_wing('M6', 10, 0.7, 5, 5)
        write_wing('NACA0012', 10, 0.7, 3, 5)
        write_wing('NACA0012', 10, 0.7, 5, 5)
