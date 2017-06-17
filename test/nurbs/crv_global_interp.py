import unittest
import numpy as np
from src.nurbs.curve import GlobalInterpolatedCrv
from src.iges.iges_core import IGES_Model


def build_airfoil(foil, p, knot_method):
    """
    构建翼型插值曲线
    :param foil: 翼型名称
    :param p: 曲线次数
    :param knot_method: 计算节点的方法
    :return: None
    """

    pnt_list = []

    fin = open('../../airfoil/' + foil + '.dat')
    model_file = IGES_Model('{}_{}_{}.igs'.format(foil, p, knot_method))

    for pnt in fin:
        x, y, z = pnt.split()
        pnt_list.append([float(x), float(y), float(z)])
    fin.close()

    pts = np.copy(pnt_list)
    foil_crv = GlobalInterpolatedCrv(pts, p, knot_method)
    model_file.add_entity(foil_crv.to_iges(1, 0, [0, 0, 1]))
    model_file.write()


class CurveGlobalInterpolationTest(unittest.TestCase):
    @staticmethod
    def test_airfoil():
        airfoil_list = ['M6', 'NACA0012', 'RAE2822']
        order_list = [3, 5]
        method_list = ['chord', 'centripetal']
        for foil in airfoil_list:
            for p in order_list:
                for method in method_list:
                    build_airfoil(foil, p, method)
