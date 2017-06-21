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


if __name__ == '__main__':
    p3d_msh = PLOT3D()
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

    LineEnding = np.array([[0, 2], [2, 0], [1, 0], [3, 1],
                           [8, 10], [10, 8], [9, 8], [11, 9],
                           [3, 11], [2, 10], [2, 6], [0, 4],
                           [1, 5], [3, 7], [10, 14], [8, 12],
                           [9, 13], [11, 15], [4, 6], [5, 4],
                           [7, 5], [12, 14], [13, 12], [15, 13],
                           [6, 14], [4, 12], [5, 13], [7, 15]])
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

    blk0_tfi_grid = Linear_TFI_3D(lambda v, w: S1(v, w),
                                  lambda v, w: sf(1.0 - v, w),
                                  lambda w, u: S3(w, u),
                                  lambda w, u: S4(w, u),
                                  lambda u, v: S6(u, v),
                                  lambda u, v: S6(u, v))

    U, V, W = 18, 26, 12
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    w_list = np.linspace(0, 1.0, W + 1)
    ppu, ppv, ppw = np.meshgrid(u_list, v_list, w_list, indexing='ij')
    blk0_tfi_grid.calc_grid(ppu, ppv, ppw)
    p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk0_tfi_grid.get_grid()))

    SS1 = Linear_TFI_2D(L[LineMap['L02']], C08, L[LineMap['L810']], L[LineMap['L210']])
    SS2 = Linear_TFI_2D(L[LineMap['L46']], L[LineMap['L412']], L[LineMap['L1214']], L[LineMap['L614']])
    SS3 = Linear_TFI_2D(C08, L[LineMap['L04']], L[LineMap['L412']], L[LineMap['L812']])
    SS4 = Linear_TFI_2D(L[LineMap['L210']], L[LineMap['L26']], L[LineMap['L614']], L[LineMap['L1014']])
    SS5 = Linear_TFI_2D(L[LineMap['L04']], L[LineMap['L02']], L[LineMap['L26']], L[LineMap['L46']])
    SS6 = Linear_TFI_2D(L[LineMap['L812']], L[LineMap['L810']], L[LineMap['L1014']], L[LineMap['L1214']])

    blk1_tfi_grid = Linear_TFI_3D(lambda v, w: SS1(v, w),
                                  lambda v, w: SS2(v, w),
                                  lambda w, u: SS3(w, u),
                                  lambda w, u: SS4(w, u),
                                  lambda u, v: SS6(u, v),
                                  lambda u, v: SS6(u, v))

    N = 15
    n_list = np.linspace(0.0, 1.0, N + 1)
    ppn, ppu, ppw = np.meshgrid(n_list, u_list, w_list, indexing='ij')
    blk1_tfi_grid.calc_grid(ppn, ppu, ppw)
    p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk1_tfi_grid.get_grid()))

    SSS1 = Linear_TFI_2D(L[LineMap['L31']], L[LineMap['L311']], L[LineMap['L119']], C19)
    SSS2 = Linear_TFI_2D(L[LineMap['L75']], L[LineMap['L715']], L[LineMap['L1513']], L[LineMap['L513']])
    SSS3 = Linear_TFI_2D(L[LineMap['L311']], L[LineMap['L37']], L[LineMap['L715']], L[LineMap['L1115']])
    SSS4 = Linear_TFI_2D(C19, L[LineMap['L15']], L[LineMap['L513']], L[LineMap['L913']])
    SSS5 = Linear_TFI_2D(L[LineMap['L37']], L[LineMap['L31']], L[LineMap['L15']], L[LineMap['L75']])
    SSS6 = Linear_TFI_2D(L[LineMap['L1115']], L[LineMap['L119']], L[LineMap['L913']], L[LineMap['L1513']])

    blk2_tfi_grid = Linear_TFI_3D(lambda v, w: SSS1(v, w),
                                  lambda v, w: SSS2(v, w),
                                  lambda w, u: SSS3(w, u),
                                  lambda w, u: SSS4(w, u),
                                  lambda u, v: SSS6(u, v),
                                  lambda u, v: SSS6(u, v))

    ppn, ppu, ppw = np.meshgrid(n_list, u_list, w_list, indexing='ij')
    blk2_tfi_grid.calc_grid(ppn, ppu, ppw)
    p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk2_tfi_grid.get_grid()))

    p3d_msh.write('{}_{}blks.xyz'.format(AirfoilName, 3))
