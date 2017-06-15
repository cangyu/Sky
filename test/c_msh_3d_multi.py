import numpy as np
import unittest
import math
from src.aircraft.frame import BWBFrame
from src.aircraft.wing import Wing
from src.nurbs.curve import Line
from src.iges.iges_core import IGES_Model


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


def bwb_geom(n, airfoil, frame_param):
    """
    构建BWB外形曲面
    :param n: 剖面数量
    :param airfoil: 剖面翼型名称
    :param frame_param: 总体描述参数
    :return: Wing Surface
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
    wg.build_sketch(3, 3)

    return wg.geom


if __name__ == '__main__':
    SectionNum = 10
    AirfoilName = 'M6'
    FrameParam = [100, 60, 20, 30, 105, 0, 45, 30]
    La = 10
    Lf = 20 * La
    Lt = 40 * La
    R = 10 * La
    Span = 15 * La

    sf, C01, C89, Cft, C08, C19 = bwb_geom(SectionNum, AirfoilName, FrameParam)
    P = np.zeros((16, 3))
    P[0] = sf(0, 0)
    P[1] = sf(1, 0)
    P[8] = sf(0, 1)
    P[9] = sf(1, 1)
    P[2] = np.copy(P[0])
    P[2][1] += R
    P[3] = np.copy(P[1])
    P[3][1] -= R
    P[10] = np.copy(P[8])
    P[10][1] += R
    P[11] = np.copy(P[9])
    P[11][1] -= R
    P[4] = np.array([La + Lt, P[0][1], 0.0], float)
    P[5] = np.array([La + Lt, P[1][1], 0.0], float)
    P[6] = np.copy(P[4])
    P[6][1] += R
    P[7] = np.copy(P[5])
    P[7][1] -= R
    P[12] = np.copy(P[8])
    P[12][0] = P[4][0]
    P[13] = np.copy(P[9])
    P[13][0] = P[5][0]
    P[14] = np.copy(P[12])
    P[14][1] = P[6][1]
    P[15] = np.copy(P[13])
    P[15][1] = P[7][1]

    LineEnding = np.array([[0, 2], [1, 0], [1, 3], [8, 10], [9, 8], [9, 11], [3, 11], [2, 10], [2, 6], [0, 4],
                           [1, 5], [3, 7], [10, 14], [8, 12], [9, 13], [11, 15], [4, 6], [5, 4], [7, 5], [12, 14], [13, 12], [15, 13], [6, 14], [4, 12], [5, 13], [7, 15]])
    LineName = []
    for i in range(len(LineEnding)):
        LineName.append("L{}{}".format(LineEnding[i][0], LineEnding[i][1]))

    LineMap = {}
    for i in range(len(LineName)):
        LineMap[LineName[i]] = i

    L = []
    for ln in LineName:
        sp, ep = LineEnding[LineMap[ln]]
        L.append(Line(P[sp], P[ep]))

    frame = IGES_Model("frame.igs")
    for line in L:
        frame.add_entity(line.to_iges())

    frame.add_entity(C08.to_iges(0, 0, [0, 0, 0]))
    frame.add_entity(C19.to_iges(0, 0, [0, 0, 0]))
    frame.add_entity(C01.to_iges(0, 0, [0, 0, 0]))
    frame.add_entity(C89.to_iges(0, 0, [0, 0, 0]))
    frame.add_entity(Cft.to_iges(0, 0, [0, 0, 0]))
    frame.write()
