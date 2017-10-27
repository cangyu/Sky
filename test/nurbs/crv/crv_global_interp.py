import unittest

import numpy as np

from settings import AIRFOIL_DIR
from src.iges import IGES_Model
from src.nurbs.curve import GlobalInterpolatedCrv

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True


class CrvGlobalInterpTest(unittest.TestCase):
    airfoil_list = ['M6', 'NACA0012', 'RAE2822']
    order_list = [3, 5]
    method_list = ['chord', 'centripetal']

    @staticmethod
    def build_airfoil(foil, p, knot_method):
        """
        构建翼型插值曲线
        :param foil: 翼型名称
        :param p: 曲线次数
        :param knot_method: 计算节点的方法
        :return: None
        """

        pnt_list = []

        fin = open(AIRFOIL_DIR + '/' + foil + '.dat')
        fn = '{}_{}_{}.igs'.format(foil, p, knot_method)
        model_file = IGES_Model(fn)

        for pnt in fin:
            x, y, z = pnt.split()
            pnt_list.append([float(x), float(y), float(z)])
        fin.close()

        pts = np.copy(pnt_list)
        foil_crv = GlobalInterpolatedCrv(pts, p, knot_method)
        model_file.add_entity(foil_crv.to_iges(1, 0, [0, 0, 1]))
        model_file.write()
        if auto_view:
            view(fn)

    @staticmethod
    def test_airfoil():
        for foil in CrvGlobalInterpTest.airfoil_list:
            for p in CrvGlobalInterpTest.order_list:
                for method in CrvGlobalInterpTest.method_list:
                    CrvGlobalInterpTest.build_airfoil(foil, p, method)
