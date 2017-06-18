import unittest
import numpy as np
from src.nurbs.curve import GlobalInterpolatedCrv
import matplotlib.pyplot as plt


def plot_airfoil_curvature(airfoil, p, method, sample_num=2000):
    """
    从上表面尾部到下表面尾部描绘翼型的曲率
    :param airfoil: 翼型名称
    :param p: 插值曲线次数
    :param method: 生成节点的方法
    :param sample_num: 画图中采样点数量
    :return: None
    """

    '''Read points'''
    fin = open('../../airfoil/' + airfoil + '.dat')
    pnt_list = []
    for pnt in fin:
        x, y, z = pnt.split()
        pnt_list.append([float(x), float(y), float(z)])
    fin.close()

    '''Construct curve'''
    pts = np.copy(pnt_list)
    crv = GlobalInterpolatedCrv(pts, p, method)

    '''Calculate curvature'''
    u = np.linspace(0, 1.0, sample_num)
    kappa = np.zeros(sample_num)
    for i in range(sample_num):
        kappa[i] = crv.curvature(u[i])

    '''Plot'''
    plt.figure()
    plt.plot(u, kappa)
    plt.xlabel('Parameter along curve')
    plt.ylabel('Curvature')
    plt.title('{}_{}_{}'.format(airfoil, p, method))
    plt.show()


class CrvCurvatureTest(unittest.TestCase):
    airfoil_list = ['M6', 'RAE2822', 'NACA0012']
    order_list = [3, 5]
    knot_method = ['chord', 'centripetal']

    @classmethod
    def test_airfoil(cls):
        for airfoil in cls.airfoil_list:
            for order in cls.order_list:
                for method in cls.knot_method:
                    plot_airfoil_curvature(airfoil, order, method)
