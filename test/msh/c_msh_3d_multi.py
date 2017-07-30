import unittest
import numpy as np
import copy
from src.aircraft.wing import Wing
from src.aircraft.frame import BWBFrame, chebshev_dist
from src.nurbs.curve import Line, Arc
from src.nurbs.surface import BilinearSurf, Coons
from src.msh.tfi import Linear_TFI_3D
from src.msh.plot3d import PLOT3D, PLOT3D_Block


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
    fn = '{}_{}blks.xyz'.format(AirfoilName, 3)
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

    C32 = Arc.from_2pnt(P[3], P[2], 180, [0, 0, -1])
    C1110 = Arc.from_2pnt(P[11], P[10], 180, [0, 0, -1])
    C23 = Arc.from_2pnt(P[2], P[3], 180, [0, 0, 1])
    C1011 = Arc.from_2pnt(P[10], P[11], 180, [0, 0, 1])

    S1 = Coons(C32, C1110, Line(P[3], P[11]), Line(P[2], P[10]))
    S3 = Coons(Line(P[3], P[11]), C19, Line(P[3], P[1]), Line(P[11], P[9]))
    S4 = Coons(Line(P[2], P[10]), C08, Line(P[2], P[0]), Line(P[10], P[8]))
    S5 = Coons(Line(P[3], P[1]), Line(P[2], P[0]), C32, C10)
    S6 = Coons(Line(P[11], P[9]), Line(P[10], P[8]), C1110, C98)

    blk0_tfi_grid = Linear_TFI_3D(lambda v, w: S1(v, w),
                                  lambda v, w: sf(1.0 - v, w),
                                  lambda w, u: S3(w, u),
                                  lambda w, u: S4(w, u),
                                  lambda u, v: S5(u, v),
                                  lambda u, v: S6(u, v))

    U, V, W = 18, 26, 12
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    w_list = np.linspace(0, 1.0, W + 1)
    ppu, ppv, ppw = np.meshgrid(u_list, v_list, w_list, indexing='ij')
    blk0_tfi_grid.calc_grid(ppu, ppv, ppw)
    p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk0_tfi_grid.get_grid()))

    SS1 = Coons(Line(P[0], P[2]), Line(P[8], P[10]), C08, Line(P[2], P[10]))
    SS2 = BilinearSurf(np.array([[P[4], P[12]], [P[6], P[14]]]))
    SS3 = Coons(C08, Line(P[4], P[12]), Line(P[0], P[4]), Line(P[8], P[12]))
    SS4 = BilinearSurf(np.array([[P[2], P[6]], [P[10], P[14]]]))
    SS5 = BilinearSurf(np.array([[P[0], P[2]], [P[4], P[6]]]))
    SS6 = BilinearSurf(np.array([[P[8], P[10]], [P[12], P[14]]]))

    blk1_tfi_grid = Linear_TFI_3D(lambda v, w: SS1(v, w),
                                  lambda v, w: SS2(v, w),
                                  lambda w, u: SS3(w, u),
                                  lambda w, u: SS4(w, u),
                                  lambda u, v: SS5(u, v),
                                  lambda u, v: SS6(u, v))

    N = 15
    n_list = np.linspace(0.0, 1.0, N + 1)
    ppn, ppu, ppw = np.meshgrid(n_list, u_list, w_list, indexing='ij')
    blk1_tfi_grid.calc_grid(ppn, ppu, ppw)
    p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk1_tfi_grid.get_grid()))

    SSS1 = Coons(Line(P[3], P[1]), Line(P[11], P[9]), Line(P[3], P[11]), C19)
    SSS2 = BilinearSurf(np.array([[P[7], P[15]], [P[5], P[13]]]))
    SSS3 = BilinearSurf(np.array([[P[3], P[7]], [P[11], P[15]]]))
    SSS4 = Coons(C19, Line(P[5], P[13]), Line(P[1], P[5]), Line(P[9], P[13]))
    SSS5 = BilinearSurf(np.array([[P[3], P[1]], [P[7], P[5]]]))
    SSS6 = BilinearSurf(np.array([[P[11], P[9]], [P[15], P[13]]]))

    blk2_tfi_grid = Linear_TFI_3D(lambda v, w: SSS1(v, w),
                                  lambda v, w: SSS2(v, w),
                                  lambda w, u: SSS3(w, u),
                                  lambda w, u: SSS4(w, u),
                                  lambda u, v: SSS5(u, v),
                                  lambda u, v: SSS6(u, v))

    ppn, ppu, ppw = np.meshgrid(n_list, u_list, w_list, indexing='ij')
    blk2_tfi_grid.calc_grid(ppn, ppu, ppw)
    p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk2_tfi_grid.get_grid()))
    p3d_msh.write(fn)
