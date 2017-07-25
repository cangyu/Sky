import unittest
import numpy as np
import math
from src.aircraft.wing import Wing, BWBFrame

try:
    from src.misc.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True


def chebshev_dist(start, end, n):
    """
    生成切比雪夫点
    :param start: 起始值
    :param end: 终止值
    :param n: 采样点数量
    :return: n个点
    """

    ang = np.linspace(math.pi, 0, n)
    pr = np.zeros(n)
    for i in range(0, n):
        pr[i] = math.cos(ang[i])
        pr[i] = start + (end - start) / 2 * (pr[i] + 1)

    return pr


def write_bwb(n, airfoil, frame_param):
    """
    BWB外形参数化建模
    :param n: 剖面数量
    :param airfoil: 剖面翼型名称
    :param frame_param: 总体描述参数 
    :return: None
    """

    gf = BWBFrame(frame_param)

    airfoil_list = []
    for i in range(n):
        airfoil_list.append(airfoil)

    thkf = np.ones(n)
    z = chebshev_dist(0, gf.Bt, n)
    xf = gf.xfront(z)
    yf = gf.yfront(z)
    xt = gf.xtail(z)
    yt = gf.ytail(z)

    wg = Wing(airfoil_list, thkf, z, xf, yf, xt, yt)

    fn = "BWB_{}_{}_{}_{}.igs".format(n, airfoil, frame_param[0], frame_param[4])
    wg.write(fn, mirror=False)

    if auto_view:
        view(fn)


class BWBWingTest(unittest.TestCase):
    @staticmethod
    def test():
        write_bwb(8, 'M6', [100, 60, 20, 30, 105, 0, 45, 30])
        write_bwb(8, 'NACA0012', [100, 60, 20, 30, 105, 0, 45, 30])


if __name__ == '__main__':
    unittest.main()
