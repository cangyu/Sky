import unittest
import numpy as np
import math
from src.grid import FluentMSH, BCType, ThomasMiddlecoff2D, LinearTFI2D, LinearTFI3D
from src.grid import single_exponential, double_exponential, hyperbolic_tangent, chebshev_dist_multi
from src.nurbs import Circle, Line, Coons, BilinearSurf
from src.wing import WingProfile


def rectangular(u, v, l, w):
    """
    矩形
    :param u: U方向网格数量
    :type u: int
    :param v: V方向网格数量
    :type v: int
    :param l: 矩形长度
    :type l: float
    :param w: 矩形宽度
    :type w: float
    :return: None
    """

    msh = LinearTFI2D(lambda u: np.array([l * u, 0, 0]),
                      lambda v: np.array([0, w * v, 0]),
                      lambda u: np.array([l * u, w, 0]),
                      lambda v: np.array([l, w * v, 0]))

    u_list = np.linspace(0, 1.0, u + 1)
    v_list = np.linspace(0, 1.0, v + 1)
    msh.calc_grid(u_list, v_list)

    fluent_grid = FluentMSH.from_str2d(msh.grid)
    fluent_grid.save('rect_{}_{}.msh'.format(u, v))


def curve_rect(U, V, L, H1, H2, H3):
    """
    曲边矩形
    :param U: U方向网格数量
    :type U: int
    :param V: V方向网格数量
    :type V: int
    :param L: 矩形长度
    :type L: float
    :param H1: 控制高度1
    :type H1: float
    :param H2: 控制高度2
    :type H2: float
    :param H3: 控制高度3
    :type H3: float
    :return: None
    """

    msh = LinearTFI2D(lambda u: np.array([u * L, 4 * H3 * u * (1 - u), 0]),
                      lambda v: np.array([0, v * H1, 0]),
                      lambda u: np.array([u * L, (H1 * (1 - u * u) + H2 * u * u), 0]),
                      lambda v: np.array([L, v * H2, 0]))

    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    msh.calc_grid(u_list, v_list)

    fluent_grid = FluentMSH.from_str2d(msh.grid)
    fluent_grid.save('crv_rect_{}_{}.msh'.format(U, V))


def airfoil(U, V, foil, ends, thk, order, arc_start, arc_end, theta, nv):
    """
    半圆翼型
    :param U: U方向网格数量
    :param V: V方向网格数量
    :param foil: 翼型名称
    :param ends: 翼型拉伸前后端点
    :param thk: 翼型厚度缩放因子
    :param order: 翼型曲线次数
    :param arc_start: 远场圆弧起始点
    :param arc_end: 远场圆弧终止点
    :param theta: 远场圆弧角度
    :param nv: 远场圆弧法向量
    :return: None
    """

    u0 = WingProfile(foil, ends, thk).nurbs_rep(order)
    u1 = Circle.from_2pnt(arc_start, arc_end, theta, nv)
    v0 = Line(u0.start, u1.start)
    v1 = Line(u0.end, u1.end)

    msh = LinearTFI2D(u0, v0, u1, v1)

    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    msh.calc_grid(u_list, v_list)

    fluent_grid = FluentMSH.from_str2d(msh.grid, bc=(BCType.Wall, BCType.Outflow, BCType.VelocityInlet, BCType.Outflow))
    fn = '{}_{}_{}.msh'.format(foil, U, V)
    fluent_grid.save(fn)


def build_airfoil_msh(foil, ending, tk, La, Lt, R, foil_order, N):
    c0 = WingProfile(foil, ending, tk).nurbs_rep(foil_order)
    c1 = Circle.from_2pnt((La, R, 0), (La, -R, 0), 180, (0, 0, 1))

    yhu = c0.start[1]
    ydu = c0.end[1]
    p = np.array([[La, yhu, 0],
                  [La, ydu, 0],
                  [La, R, 0],
                  [La, -R, 0],
                  [La + Lt, yhu, 0],
                  [La + Lt, ydu, 0],
                  [La + Lt, R, 0],
                  [La + Lt, -R, 0]])

    c2 = Line(p[0], p[2])
    c3 = Line(p[1], p[3])
    c4 = Line(p[4], p[6])
    c5 = Line(p[5], p[7])
    c6 = Line(p[4], p[0])
    c7 = Line(p[6], p[2])
    c8 = Line(p[5], p[1])
    c9 = Line(p[7], p[3])
    c10 = Line(p[1], p[0])
    c11 = Line(p[5], p[4])

    pu = [double_exponential(N[0], 0.5, -1.5, 0.5),  # c0, c1
          hyperbolic_tangent(N[1], 2),  # c2, c3, c4, c5
          single_exponential(N[2], -3),  # c6, c7, c8, c9
          np.linspace(0.0, 1.0, N[3])]  # c10, c11

    '''翼型前缘'''
    grid1 = LinearTFI2D(c0, c2, c1, c3)
    grid1.calc_grid(pu[0], pu[1])
    tm_grid1 = ThomasMiddlecoff2D(grid1.get_grid())
    tm_grid1.calc_grid()

    '''翼型后缘上部'''
    grid2 = LinearTFI2D(c6, c4, c7, c2)
    grid2.calc_grid(pu[2], pu[1])

    '''翼型后缘下部'''
    grid3 = LinearTFI2D(c8, c5, c9, c3)
    grid3.calc_grid(pu[2], pu[1])

    '''钝尾缘部分'''
    grid4 = LinearTFI2D(c8, c11, c6, c10)
    grid4.calc_grid(pu[2], pu[3])

    '''网格, 边界条件, 邻接关系'''
    blk = [tm_grid1.get_grid(), grid2.get_grid(), grid3.get_grid(), grid4.get_grid()]

    bc = [(BCType.Wall, BCType.Interior, BCType.PressureFarField, BCType.Interior),
          (BCType.Interior, BCType.PressureFarField, BCType.PressureFarField, BCType.Interior),
          (BCType.Interior, BCType.PressureFarField, BCType.PressureFarField, BCType.Interior),
          (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Wall)]

    adj = [((0, 2), (1, 4)),
           ((2, 4), (0, 4)),
           ((3, 3), (1, 1)),
           ((2, 1), (3, 1)),
           ((0, 3), (0, 0)),
           ((0, 0), (0, 1)),
           ((1, 3), (0, 0)),
           ((0, 0), (2, 3)),
           ((0, 0), (3, 4)),
           ((3, 2), (0, 0)),
           ((1, 2), (0, 0)),
           ((0, 0), (2, 2))]

    '''构建MSH文件'''
    fn = '{}_C_grid.msh'.format(foil)
    msh = FluentMSH.from_str2d_multi(blk, bc, adj)
    msh.save(fn)


def sect(r_min, r_max, theta_min, theta_max, h_min, h_max, nu, nv, nw):
    """
    扇柱单块网格
    :param r_min: 内径
    :type r_min: float
    :param r_max: 外径
    :type r_max: float
    :param theta_min: 起始角
    :type theta_min: float
    :param theta_max: 终止角
    :type theta_max: float
    :param h_min: 底部高度
    :type h_min: float
    :param h_max: 顶部高度
    :type h_max: float
    :param nu: X方向网格数
    :type nu: int
    :param nv: Y方向网格数
    :type nv: int
    :param nw: Z方向网格数
    :type nw: int
    :return: TFI_3D Mesh
    """

    def share(p, a, b):
        return (1 - p) * a + p * b

    def s1(v, w):
        return np.array([r_min * math.cos(math.radians(share(v, theta_min, theta_max))),
                         r_min * math.sin(math.radians(share(v, theta_min, theta_max))),
                         share(w, h_min, h_max)])

    def s2(v, w):
        return np.array([r_max * math.cos(math.radians(share(v, theta_min, theta_max))),
                         r_max * math.sin(math.radians(share(v, theta_min, theta_max))),
                         (1 - w) * h_min + w * h_max])

    def s3(w, u):
        return np.array([share(u, r_min, r_max) * math.cos(math.radians(theta_min)),
                         share(u, r_min, r_max) * math.sin(math.radians(theta_min)),
                         share(w, h_min, h_max)])

    def s4(w, u):
        return np.array([share(u, r_min, r_max) * math.cos(math.radians(theta_max)),
                         share(u, r_min, r_max) * math.sin(math.radians(theta_max)),
                         share(w, h_min, h_max)])

    def s5(u, v):
        return np.array([share(u, r_min, r_max) * math.cos(math.radians(share(v, theta_min, theta_max))),
                         share(u, r_min, r_max) * math.sin(math.radians(share(v, theta_min, theta_max))),
                         h_min])

    def s6(u, v):
        return np.array([share(u, r_min, r_max) * math.cos(math.radians(share(v, theta_min, theta_max))),
                         share(u, r_min, r_max) * math.sin(math.radians(share(v, theta_min, theta_max))),
                         h_max])

    msh = LinearTFI3D(s1, s2, s3, s4, s5, s6)

    u_list = np.linspace(0, 1.0, nu + 1)
    v_list = np.linspace(0, 1.0, nv + 1)
    w_list = np.linspace(0, 1.0, nw + 1)
    msh.calc_grid(u_list, v_list, w_list)

    return msh.grid


class FluentMSHTestCase(unittest.TestCase):
    def test_2d_single(self):
        rectangular(30, 10, 3, 1)
        curve_rect(50, 25, 100, 40, 60, 10)
        airfoil(130, 45, 'M6', [[0, 0, 5], [20, 0, 5]], 1.6, 5, [25, 60, 5], [25, -60, 5], 180, [0, 0, 1])
        self.assertTrue(True)

    def test_2d_multi(self):
        build_airfoil_msh('M6', [[0, 0, 0], [10, 0, 0]], 1.2, 10, 200, 100, 5, [41, 26, 61, 2])
        build_airfoil_msh('NACA0012', [[0, 0, 0], [10, 0, 0]], 1.2, 10, 200, 100, 5, [41, 26, 61, 2])
        self.assertTrue(True)

    def test_3d_single(self):
        grid = sect(10, 32, 0, 160, -15, 15, 25, 30, 20)
        fluent_msh = XF_MSH.from_str3d(grid)
        fluent_msh.save("sect.msh")
        self.assertTrue(True)

    def test_3d_multi(self):
        self.assertTrue(True)

        c_root = 10
        c_mid = 6
        c_tip = 1.1
        b_mid = 4.5
        b_tip = 13
        alpha_mid = 35
        alpha_tip = 30
        frm = BWBFrame(c_root, c_mid, c_tip, b_mid, b_tip, alpha_mid, alpha_tip)

        inner_sec_num = 6
        outer_sec_num = 8
        u_mid = b_mid / b_tip
        u_dist = chebshev_dist_multi([0, u_mid, 1], [inner_sec_num, outer_sec_num])
        sec_num = len(u_dist)

        foil = ['NACA0012', 'NACA0012', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6']
        z_offset = list(map(frm.z, u_dist))
        length = list(map(lambda u: frm.x_tail(u) - frm.x_front(u), u_dist))
        sweep_back = list(map(lambda u: math.degrees(math.atan2(frm.x_front(u), frm.z(u))), u_dist))
        twist = np.array([3, 2.5, 2.2, 2, 2, 1.8, 1.6, 1.1, 1, 1, 0.5, 0.2, 0])
        dihedral = np.array([0, 1, 1.1, 1.2, 2, 2, 2, 2, 2, 2, 2.5, 3, 3])
        twist_pos = np.full(sec_num, 0.25)
        y_ref = np.zeros(sec_num)
        thickness_factor = np.ones(sec_num)

        wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)

        La = length[0]
        Lt = 40 * La
        R = 20 * La
        Span = z_offset[-1]

        sf = wg.surf()
        C01 = Line(wg.section[0].tail_up, wg.section[0].tail_down)
        C89 = Line(wg.section[-1].tail_up, wg.section[-1].tail_down)
        C08 = sf.extract('U', 0)
        C19 = sf.extract('U', 1)

        C10 = deepcopy(C01)
        C10.reverse()
        C98 = deepcopy(C89)
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

        C32 = Circle.from_2pnt(P[3], P[2], 180, [0, 0, -1])
        C1110 = Circle.from_2pnt(P[11], P[10], 180, [0, 0, -1])
        C23 = Circle.from_2pnt(P[2], P[3], 180, [0, 0, 1])
        C1011 = Circle.from_2pnt(P[10], P[11], 180, [0, 0, 1])

        S1 = Coons(C23, C1011, Line(P[2], P[10]), Line(P[3], P[11]))
        S3 = Coons(C08, Line(P[2], P[10]), Line(P[0], P[2]), Line(P[8], P[10]))
        S4 = Coons(C19, Line(P[3], P[11]), Line(P[1], P[3]), Line(P[9], P[11]))
        S5 = Coons(Line(P[0], P[2]), Line(P[1], P[3]), C01, C23)
        S6 = Coons(Line(P[8], P[10]), Line(P[9], P[11]), C89, C1011)

        blk0_tfi_grid = LinearTFI3D(lambda v, w: sf(v, w),
                                    lambda v, w: S1(v, w),
                                    lambda w, u: S3(w, u),
                                    lambda w, u: S4(w, u),
                                    lambda u, v: S5(u, v),
                                    lambda u, v: S6(u, v))

        U, V, W = 54, 62, 40
        u_list = np.linspace(0, 1.0, U + 1)
        v_list = np.linspace(0, 1.0, V + 1)
        w_list = np.linspace(0, 1.0, W + 1)
        blk0_tfi_grid.calc_grid(u_list, v_list, w_list)

        SS1 = Coons(Line(P[0], P[2]), Line(P[8], P[10]), C08, Line(P[2], P[10]))
        SS2 = BilinearSurf(np.array([[P[4], P[12]], [P[6], P[14]]]))
        SS3 = Coons(C08, Line(P[4], P[12]), Line(P[0], P[4]), Line(P[8], P[12]))
        SS4 = BilinearSurf(np.array([[P[2], P[6]], [P[10], P[14]]]))
        SS5 = BilinearSurf(np.array([[P[0], P[2]], [P[4], P[6]]]))
        SS6 = BilinearSurf(np.array([[P[8], P[10]], [P[12], P[14]]]))

        blk1_tfi_grid = LinearTFI3D(lambda v, w: SS1(v, w),
                                    lambda v, w: SS2(v, w),
                                    lambda w, u: SS3(w, u),
                                    lambda w, u: SS4(w, u),
                                    lambda u, v: SS5(u, v),
                                    lambda u, v: SS6(u, v))

        N = 60
        n_list = np.linspace(0.0, 1.0, N + 1)
        blk1_tfi_grid.calc_grid(n_list, u_list, w_list)

        SSS1 = Coons(Line(P[1], P[5]), Line(P[9], P[13]), C19, Line(P[5], P[13]))
        SSS2 = BilinearSurf(np.array([[P[3], P[11]], [P[7], P[15]]]))
        SSS3 = Coons(C19, Line(P[3], P[11]), Line(P[1], P[3]), Line(P[9], P[11]))
        SSS4 = BilinearSurf(np.array([[P[5], P[7]], [P[13], P[15]]]))
        SSS5 = BilinearSurf(np.array([[P[1], P[5]], [P[3], P[7]]]))
        SSS6 = BilinearSurf(np.array([[P[9], P[13]], [P[11], P[15]]]))

        blk2_tfi_grid = LinearTFI3D(lambda v, w: SSS1(v, w),
                                    lambda v, w: SSS2(v, w),
                                    lambda w, u: SSS3(w, u),
                                    lambda w, u: SSS4(w, u),
                                    lambda u, v: SSS5(u, v),
                                    lambda u, v: SSS6(u, v))

        blk2_tfi_grid.calc_grid(u_list, n_list, w_list)

        SSSI1 = Coons(Line(P[0], P[4]), Line(P[8], P[12]), C08, Line(P[4], P[12]))
        SSSI2 = Coons(Line(P[1], P[5]), Line(P[9], P[13]), C19, Line(P[5], P[13]))
        SSSI3 = Coons(C08, C19, Line(P[0], P[1]), Line(P[8], P[9]))
        SSSI4 = BilinearSurf(np.array([[P[4], P[5]], [P[12], P[13]]]))
        SSSI5 = BilinearSurf(np.array([[P[0], P[4]], [P[1], P[5]]]))
        SSSI6 = BilinearSurf(np.array([[P[8], P[12]], [P[9], P[13]]]))

        blk3_tfi_grid = LinearTFI3D(lambda v, w: SSSI1(v, w),
                                    lambda v, w: SSSI2(v, w),
                                    lambda w, u: SSSI3(w, u),
                                    lambda w, u: SSSI4(w, u),
                                    lambda u, v: SSSI5(u, v),
                                    lambda u, v: SSSI6(u, v))

        M = 16
        m_list = np.linspace(0, 1, M + 1)
        blk3_tfi_grid.calc_grid(m_list, n_list, w_list)

        '''网格, 边界条件, 邻接关系'''
        blk = [blk0_tfi_grid.get_grid(), blk1_tfi_grid.get_grid(), blk2_tfi_grid.get_grid(), blk3_tfi_grid.get_grid()]

        bc = [(BCType.Wall, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Symmetry, BCType.Wall),
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Wall),
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Wall),
              (BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField, BCType.Symmetry, BCType.Wall)]

        adj = [((0, 0), (0, 1), 0, False),
               ((0, 2), (0, 0), 0, False),
               ((1, 1), (0, 3), 1, True),
               ((0, 4), (2, 3), 0, False),
               ((0, 0), (0, 5), 0, False),
               ((0, 6), (0, 0), 0, False),

               ((1, 2), (0, 0), 0, False),
               ((3, 1), (1, 3), 1, True),
               ((1, 4), (0, 0), 0, False),
               ((0, 0), (1, 5), 0, False),
               ((1, 6), (0, 0), 0, False),

               ((3, 2), (2, 1), 1, False),
               ((2, 2), (0, 0), 0, False),
               ((2, 4), (0, 0), 0, False),
               ((0, 0), (2, 5), 0, False),
               ((2, 6), (0, 0), 0, False),

               ((0, 0), (3, 3), 0, False),
               ((3, 4), (0, 0), 0, False),
               ((0, 0), (3, 5), 0, False),
               ((3, 6), (0, 0), 0, False)]

        '''构建MSH文件'''
        msh = XF_MSH.from_str3d_multi(blk, bc, adj)
        msh.save('3d_multi_blk.msh')


if __name__ == '__main__':
    unittest.main()
