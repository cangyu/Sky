import numpy as np
import unittest
import math
import copy
from src.aircraft.frame import BWBFrame
from src.aircraft.wing import Wing
from src.nurbs.curve import Line, Arc
from src.iges.iges_core import IGES_Model
from src.com.catia import view
from src.msh.tfi import Linear_TFI_2D, Linear_TFI_3D
from src.msh.plot3d import PLOT3D, PLOT3D_Block


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


def write_uniform_p3d(msh: Linear_TFI_3D, U: int, V: int, W: int, fn="msh_p3d.xyz"):
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    w_list = np.linspace(0, 1.0, W + 1)
    ppu, ppv, ppw = np.meshgrid(u_list, v_list, w_list, indexing='ij')
    msh.calc_grid(ppu, ppv, ppw)
    grid = PLOT3D()
    grid.add_block(PLOT3D_Block.build_from_3d(msh.get_grid()))
    grid.write(fn)


if __name__ == '__main__':
    SectionNum = 10
    AirfoilName = 'M6'
    FrameParam = [100, 60, 20, 30, 105, 0, 45, 30]
    La = 10
    Lt = 40 * La
    R = 20 * La
    Span = 15 * La

    sf, C01, C89, Cft, C08, C19 = bwb_geom(SectionNum, AirfoilName, FrameParam)
    C10 = copy.deepcopy(C01)
    C10.reverse()
    C98 = copy.deepcopy(C89)
    C98.reverse()

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

    LineEnding = np.array([[2, 0], [1, 0], [3, 1], [10, 8], [9, 8], [11, 9], [3, 11], [2, 10], [2, 6], [0, 4],
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

    C32 = Arc.from_2pnt(P[3], P[2], 180, [0, 0, -1])
    C1110 = Arc.from_2pnt(P[11], P[10], 180, [0, 0, -1])
    C23 = Arc.from_2pnt(P[2], P[3], 180, [0, 0, 1])
    C1011 = Arc.from_2pnt(P[10], P[11], 180, [0, 0, 1])

    frame = IGES_Model("frame.igs")
    for line in L:
        frame.add_entity(line.to_iges())

    frame.add_entity(C08.to_iges(0, 0, [0, 0, 0]))
    frame.add_entity(C19.to_iges(0, 0, [0, 0, 0]))
    frame.add_entity(C01.to_iges(0, 0, [0, 0, 0]))
    frame.add_entity(C89.to_iges(0, 0, [0, 0, 0]))
    frame.add_entity(Cft.to_iges(0, 0, [0, 0, 0]))
    frame.add_entity(C23.to_iges(1, 0, [0, 0, 1]))
    frame.add_entity(C1011.to_iges(1, 0, [0, 0, 1]))
    frame.write()

    # view(frame.filename)

    S1 = Linear_TFI_2D(C32, L[LineMap['L311']], C1110, L[LineMap['L210']])
    S3 = Linear_TFI_2D(L[LineMap['L311']], L[LineMap['L31']], C19, L[LineMap['L119']])
    S4 = Linear_TFI_2D(L[LineMap['L210']], L[LineMap['L20']], C08, L[LineMap['L108']])
    S5 = Linear_TFI_2D(L[LineMap['L31']], C32, L[LineMap['L20']], C10)
    S6 = Linear_TFI_2D(L[LineMap['L119']], C1110, L[LineMap['L108']], C98)

    tfi_grid = Linear_TFI_3D(lambda v, w: S1(v, w),
                             lambda v, w: sf(1.0 - v, w),
                             lambda w, u: S3(w, u),
                             lambda w, u: S4(w, u),
                             lambda u, v: S6(u, v),
                             lambda u, v: S6(u, v))

    write_uniform_p3d(tfi_grid, 60, 140, 50)
