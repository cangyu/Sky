import unittest
import os
import math
import numpy as np
from numpy.linalg import norm
from copy import deepcopy
import matplotlib.pyplot as plt
from tfi import LinearTFI2D, LinearTFI3D
from plot3d import Plot3D, Plot3DBlock
from smooth import Laplace2D, ThomasMiddlecoff2D
from fluent import XF_MSH, BCType
from misc import pnt_dist, read_airfoil_pts
from spacing import hyperbolic_tangent, uniform, single_exponential, double_exponential
from iges import Model, Entity116
from nurbs import Crv, Circle, Line, Spline, ConicArc, Surf, Skinned, RuledSurf
from settings import AIRFOIL_LIST


class Airfoil(object):
    def __init__(self, foil):
        """
        2D Airfoil, with chord length equals to 1.
        :param foil: Airfoil name(Capital Case).
        :type foil: str
        """

        if foil not in AIRFOIL_LIST:
            raise FileNotFoundError('Airfoil \'{}\' not included at present.'.format(foil))

        self.name = foil
        self.pts = read_airfoil_pts(foil)

    def __repr__(self):
        return '{} with {} points'.format(self.name, self.pnt_num)

    @property
    def z(self):
        return self.pts[0][2]

    @property
    def pnt_num(self):
        return len(self.pts)

    @property
    def crv(self):
        return Spline(self.pts)

    @property
    def tail_up(self):
        return self.pts[0]

    @property
    def tail_down(self):
        return self.pts[-1]

    @property
    def tail(self):
        return (self.tail_up + self.tail_down) / 2

    @property
    def front(self):
        """
        The most front point of the airfoil.
        """

        total = self.pnt_num
        cx = self.pts[0][0]
        k = 1
        while k < total and self.pts[k][0] < cx:
            cx = self.pts[k][0]
            k += 1

        return self.pts[k - 1]

    @property
    def chord_len(self):
        return pnt_dist(self.front, self.tail)

    @property
    def is_blunt(self):
        return not math.isclose(norm(self.tail_up - self.tail_down), 0)

    def to_blunt(self):
        pfx = self.front[0]
        r0 = self.chord_len
        self.pts = self.pts[1:-1]
        r1 = self.chord_len
        ratio = r0 / r1
        for k in range(len(self.pts)):
            self.pts[k][0] = pfx + (self.pts[k][0] - pfx) * ratio

    def save(self, fn):
        """
        Save all coordinates into file.
        :param fn: File name.
        :type fn: str
        :return: None.
        """

        f_out = open(fn, 'w')
        for p in self.pts:
            f_out.write('{:10.6f}\t{:10.6f}\t{:10.6f}\n'.format(p[0], p[1], p[2]))
        f_out.close()

    def show(self):
        (px, py, pz) = zip(*self.pts)
        plt.plot(px, py)
        plt.gca().set_aspect('equal')
        plt.show()

    def plot_curvature(self, n=2000):
        """
        Plot the curvature along the airfoil from tail-up to tail_down, anti-clockwise.
        By default, the interpolated curve is cubic
        :param n: Number of sampling points.
        :type n: int
        :return: None
        """

        '''Build curve and calculate curvature'''
        crv = self.crv
        ul = np.linspace(crv.U[0], crv.U[-1], n)
        kappa = list(map(lambda u: crv.curvature(u), ul))

        '''Plot'''
        plt.figure()
        plt.plot(ul, kappa)
        plt.xlabel('Parameter along curve')
        plt.ylabel('Curvature')
        plt.title(self.name)
        plt.show()

    def gen_grid(self, *args, **kwargs):
        """
        Generate grid for 2D airfoil or wing profile.
        :param args: Containing the geometric description and node distribution of the grid.
        :param kwargs: Extra options on smooth and spacing.
        :return: The wire-frame, plot3d-grid and fluent-grid(with predefined BC) of the flow field.
        """

        '''
        a: Width of the front part of the flow field.
        b: Semi-height of the flow field.
        c: Width of the rear part of the flow field.
        n0: Num of points along airfoil.
        n1: Num of points along vertical direction.
        n2: Num of points along horizontal direction in rear field.
        n3: Num of Points on the trailing edge.
        '''
        assert len(args) == 7
        a, b, c, n0, n1, n2, n3 = args

        wire_frame = Model()
        p3d_grid = Plot3D()

        '''Flow-field geometries'''
        pts = np.empty((8, 3), float)
        pts[0] = self.tail_up
        pts[1] = self.tail_down
        pts[2] = np.array([0, b, self.z])
        pts[3] = np.array([0, -b, self.z])
        pts[4] = np.array([c, pts[0][1], self.z])
        pts[5] = np.array([c, pts[1][1], self.z])
        pts[6] = np.array([c, b, self.z])
        pts[7] = np.array([c, -b, self.z])

        crv = [self.crv,  # c0
               ConicArc(pts[2], (-1, 0, 0), pts[3], (1, 0, 0), (-a, 0, 0)),  # c1
               Line(pts[0], pts[2]),  # c2
               Line(pts[1], pts[3]),  # c3
               Line(pts[4], pts[6]),  # c4
               Line(pts[5], pts[7]),  # c5
               Line(pts[0], pts[4]),  # c6
               Line(pts[1], pts[5]),  # c7
               Line(pts[2], pts[6]),  # c8
               Line(pts[3], pts[7]),  # c9
               Line(pts[0], pts[1]),  # c10
               Line(pts[4], pts[5])]  # c11

        '''Construct wire-frame'''
        for p in pts:
            wire_frame.add(Entity116(p[0], p[1], p[2]))
        for c in crv:
            wire_frame.add(c.to_iges())

        '''Knot distribution'''
        u = [double_exponential(n0, 0.5, -1.5, 0.5),  # c0, c1
             hyperbolic_tangent(n1, 2),  # c2, c3, c4, c5
             single_exponential(n2, 3),  # c6, c7, c8, c9
             uniform(n3)]  # c10, c11

        '''Structured grid blocks'''
        leading_blk = LinearTFI2D(crv[2], crv[0], crv[3], crv[1])
        tailing_up_blk = LinearTFI2D(crv[6], crv[2], crv[8], crv[4])
        tailing_down_blk = LinearTFI2D(crv[3], crv[7], crv[5], crv[9])
        rear_blk = LinearTFI2D(crv[10], crv[6], crv[11], crv[7])

        '''Construct Plot3D grid for basic checking'''
        leading_blk.calc_grid(u[1], u[0])
        leading_tfi_grid = leading_blk.grid
        leading_smooth_ok = False
        if 'leading_smooth' in kwargs:
            smooth = kwargs['leading_smooth']
            if smooth in ('Laplace', 'laplace'):
                leading_grid_laplace = Laplace2D(leading_tfi_grid)
                leading_grid_laplace.smooth()
                p3d_grid.add(Plot3DBlock.construct_from_array(leading_grid_laplace.grid))
                leading_smooth_ok = True
            if smooth in ('TM', 'tm', 'Thomas-Middlecoff', 'thomas-middlecoff'):
                leading_grid_tm = ThomasMiddlecoff2D(leading_tfi_grid)
                leading_grid_tm.smooth()
                p3d_grid.add(Plot3DBlock.construct_from_array(leading_grid_tm.grid))
                leading_smooth_ok = True
        if not leading_smooth_ok:
            p3d_grid.add(Plot3DBlock.construct_from_array(leading_tfi_grid))

        tailing_up_blk.calc_grid(u[2], u[1])
        tailing_up_grid = tailing_up_blk.grid
        p3d_grid.add(Plot3DBlock.construct_from_array(tailing_up_grid))

        tailing_down_blk.calc_grid(u[1], u[2])
        tailing_down_grid = tailing_down_blk.grid
        p3d_grid.add(Plot3DBlock.construct_from_array(tailing_down_grid))

        if self.is_blunt:
            rear_blk.calc_grid(u[3], u[2])
            rear_grid = rear_blk.grid
            p3d_grid.add(Plot3DBlock.construct_from_array(rear_grid))

        return wire_frame, p3d_grid


class WingProfile(Airfoil):
    INTRINSIC_PARAM = ['Airfoil', 'Z(m)', 'X_front(m)', 'Y_front(m)', 'X_tail(m)', 'Y_tail(m)', 'Thickness Ratio']
    GEOM_PARAM = ['Airfoil', 'Z(m)', 'Length(m)', 'SweepBack(deg)', 'Twist(deg)', 'Dihedral(deg)', 'TwistPos', 'Thickness Ratio']

    def __init__(self, foil, ends, thickness_factor=1.0):
        """
        3D profile at certain position.
        :param foil: Airfoil name.
        :type foil: str
        :param ends: Starting and ending points of the profile.
        :param thickness_factor: Vertical stretching factor.
        :type thickness_factor: float
        """

        super(WingProfile, self).__init__(foil)

        '''Inspect endings'''
        self.ending = np.copy(ends)
        if not math.isclose(ends[0][2], ends[1][2]):
            raise AssertionError("Invalid ending coordinates in Z direction!")

        cl = self.chord_len
        if math.isclose(cl, 0.0):
            raise ZeroDivisionError("Invalid ending coordinates in XY direction!")

        rotation = complex((ends[1][0] - ends[0][0]) / cl, (ends[1][1] - ends[0][1]) / cl)

        '''Build profile'''
        if not self.is_blunt:
            self.to_blunt()
        for i in range(self.pnt_num):
            '''Stretch, Z offset and Thickness'''
            self.pts[i][0] *= cl
            self.pts[i][1] *= (cl * thickness_factor)
            self.pts[i][2] = ends[0][2]

            '''Rotate around ends[0]'''
            origin_vector = complex(self.pts[i][0], self.pts[i][1])
            origin_vector *= rotation
            self.pts[i][0] = origin_vector.real
            self.pts[i][1] = origin_vector.imag

            '''Move to ends[0]'''
            self.pts[i][0] += ends[0][0]
            self.pts[i][1] += ends[0][1]

    @property
    def front(self):
        return self.ending[0]

    @property
    def tail(self):
        return self.ending[-1]

    @classmethod
    def from_geom_param(cls, foil, z_offset, length, sweep_back, twist, dihedral, twist_pos=0.25, y_ref=0, thickness_factor=1.0):
        """
        从几何描述参数构建机翼剖面
        :param foil: 翼型名称
        :type foil: str
        :param z_offset: Z方向偏移量
        :type z_offset: float
        :param length: 剖面长度
        :type length: float
        :param sweep_back: 后掠角
        :type sweep_back: float
        :param twist: 相对翼根弦线的扭转角(默认在1/4弦长处扭转)
        :type twist: float
        :param dihedral: 相对翼根的上反角
        :type dihedral: float
        :param twist_pos: 扭转中心
        :type twist_pos: float
        :param y_ref: 翼根处Y方向基准坐标
        :type y_ref: float
        :param thickness_factor: 纵向厚度拉伸系数
        :type thickness_factor: float
        :return: 机翼剖面
        :rtype: WingProfile
        """

        x_offset = z_offset * math.tan(math.radians(sweep_back))
        y_offset = y_ref + z_offset * math.tan(math.radians(dihedral))
        front = np.array([x_offset, y_offset, z_offset], float)
        tail = np.array([x_offset + length, y_offset, z_offset], float)
        center = (1 - twist_pos) * front + twist_pos * tail
        theta = math.radians(-twist)
        rot = complex(math.cos(theta), math.sin(theta))
        d1 = front - center
        d2 = tail - center
        c1 = complex(d1[0], d1[1]) * rot
        c2 = complex(d2[0], d2[1]) * rot
        front[0] = center[0] + c1.real
        front[1] = center[1] + c1.imag
        tail[0] = center[0] + c2.real
        tail[1] = center[1] + c2.imag
        ending = np.array([front, tail])
        return cls(foil, ending, thickness_factor)


class Wing(object):
    def __init__(self, profiles):
        """
        从剖面序列构造机翼
        :param profiles: 机翼剖面序列
        """

        self.profile = deepcopy(profiles)

    def __repr__(self):
        return "Wing with {} sections".format(self.size)

    @property
    def size(self):
        """
        Number of profiles within this wing.
        :return: Num of profiles.
        :rtype: int
        """

        return len(self.profile)

    @property
    def root(self):
        """
        Profile at root.
        :return: The root profile in WingProfile representation.
        :rtype: WingProfile
        """

        return self.profile[0]

    @property
    def tip(self):
        """
        Profile at tip.
        :return: The tip profile in WingProfile representation.
        :rtype: WingProfile
        """

        return self.profile[-1]

    @property
    def surf(self):
        """
        构建机翼轮廓曲线、曲面
        :return: 机翼蒙皮曲面
        :rtype: Skinned
        """

        profile_list = [self.profile[i].crv for i in range(self.size)]
        return Skinned(profile_list, 3, 3)

    @property
    def leading(self):
        """
        Leading edge of the wing.
        :return: Leading edge in NURBS representation.
        :rtype: Crv
        """

        pts = [self.profile[i].front for i in range(self.size)]
        return Spline(pts, method='chord')

    @property
    def tailing_up(self):
        return self.surf.extract('U', 0)

    @property
    def tailing_down(self):
        return self.surf.extract('U', 1)

    def iges_model(self, mirror=True):
        """
        生成机翼相应的IGES模型
        :param mirror: 是否生成对称部分
        :type mirror: bool
        :return: 可用于后续生成IGS文件的IGES_Model对象
        :rtype: Model
        """

        wing_model = Model()

        '''Leading edge'''
        wing_model.add(self.leading.to_iges())

        '''Profiles'''
        for elem in self.profile:
            wing_model.add(elem.crv.to_iges())

        '''Skin'''
        sk = self.surf
        wing_model.add(sk.to_iges())

        '''Surf mirror'''
        if mirror:
            msk = deepcopy(sk)
            for i in range(msk.n + 1):
                for j in range(msk.m + 1):
                    msk.Pw[i][j][2] *= -1
            wing_model.add(msk.to_iges())

        return wing_model

    @classmethod
    def from_intrinsic_desc(cls, airfoil, thickness, z, xf, yf, xt, yt):
        n = len(airfoil)
        section_list = []
        for k in range(n):
            ends = np.empty((2, 3), float)
            ends[0][2] = ends[1][2] = z[k]
            ends[0][0] = xf[k]
            ends[0][1] = yf[k]
            ends[1][0] = xt[k]
            ends[1][1] = yt[k]
            section_list.append(WingProfile(airfoil[k], ends, thickness[k]))

        return cls(section_list)

    @classmethod
    def from_geom_desc(cls, airfoil, length, thickness, z, sweep, twist, twist_pos, dihedral, y_ref):
        n = len(airfoil)
        section_list = []
        for k in range(n):
            section_list.append(WingProfile.from_geom_param(airfoil[k], z[k], length[k], sweep[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness[k]))

        return cls(section_list)

    @classmethod
    def from_frame(cls, airfoil, thickness, u, frm):
        """
        根据给定的参数化模型生成机翼
        :param airfoil: 剖面翼型序列
        :param thickness: 剖面厚度拉伸系数
        :param u: 剖面位置分布参数
        :param frm: 参数化模型
        :type frm: WingFrame
        :return:
        """

        z = list(map(frm.z, u))
        xf = list(map(frm.x_front, u))
        yf = list(map(frm.y_front, u))
        xt = list(map(frm.x_tail, u))
        yt = list(map(frm.y_tail, u))
        return cls.from_intrinsic_desc(airfoil, thickness, z, xf, yf, xt, yt)

    def gen_grid(self, fn, enn=50, inner_brk=(0.44, 0.56), outer_brk=(0.3, 0.72)):
        """
        Generate the multi-block grid for a simplified wing.
        :return: None
        """

        '''Grid parameters'''
        la = self.profile[0].chord_len
        lt = 30 * la
        r = 10 * la
        inner_spn = self.profile[-1].z_offset
        outer_spn = 20 * inner_spn

        wsf = self.surf()
        crv_root = wsf.extract('V', 0)
        crv_tip = wsf.extract('V', 1)
        far = WingProfile.from_geom_param(self.profile[-1].name, inner_spn + outer_spn, la, 0, 0, 0, thickness_factor=3)
        crv_far = far.crv
        fsf = RuledSurf(crv_tip, crv_far)

        brk_root = brk_tip = brk_far = inner_brk
        obrk_root = obrk_tip = obrk_far = outer_brk

        '''Points, lines, curves, surfs'''
        p = np.zeros((36, 3))
        p[2] = wsf(0, 0)
        p[4] = wsf(1, 0)
        p[0] = p[2]
        p[0][1] += r
        p[6] = p[4]
        p[6][1] -= r
        p[1] = p[0]
        p[1][0] += lt
        p[3] = p[2]
        p[3][0] += lt
        p[5] = p[4]
        p[5][0] += lt
        p[7] = p[6]
        p[7][0] += lt
        p[9] = p[1]
        p[9][2] += inner_spn
        p[11] = p[3]
        p[11][2] += inner_spn
        p[13] = p[5]
        p[13][2] += inner_spn
        p[15] = p[7]
        p[15][2] += inner_spn
        p[8] = p[0]
        p[8][2] += inner_spn
        p[10] = wsf(0, 1)
        p[12] = wsf(1, 1)
        p[14] = p[6]
        p[14][2] += inner_spn
        p[16] = p[8]
        p[16][2] += outer_spn
        p[17] = p[9]
        p[17][2] += outer_spn
        p[22] = p[14]
        p[22][2] += outer_spn
        p[23] = p[15]
        p[23][2] += outer_spn
        p[18] = crv_far.start
        p[20] = crv_far.end
        p[19] = p[18]
        p[19][0] += lt
        p[21] = p[20]
        p[21][0] += lt
        p[24] = crv_root(brk_root[0])
        p[25] = crv_root(brk_root[1])
        p[26] = crv_tip(brk_tip[0])
        p[27] = crv_tip(brk_tip[1])
        p[28] = crv_far(brk_far[0])
        p[29] = crv_far(brk_far[1])
        outer_root = Circle.from_2pnt(p[0], p[6], 180, (0, 0, 1))
        outer_tip = Circle.from_2pnt(p[8], p[14], 180, (0, 0, 1))
        outer_far = Circle.from_2pnt(p[16], p[22], 180, (0, 0, 1))
        p[30] = outer_root(obrk_root[0])
        p[31] = outer_root(obrk_root[1])
        p[32] = outer_tip(obrk_tip[0])
        p[33] = outer_tip(obrk_tip[1])
        p[34] = outer_far(obrk_far[0])
        p[35] = outer_far(obrk_far[1])

        l = [Line(p[0], p[1]),  # 0
             Line(p[2], p[3]),  # 1
             Line(p[4], p[5]),  # 2
             Line(p[6], p[7]),  # 3
             Line(p[8], p[9]),  # 4
             Line(p[10], p[11]),  # 5
             Line(p[12], p[13]),  # 6
             Line(p[14], p[15]),  # 7
             Line(p[16], p[17]),  # 8
             Line(p[18], p[19]),  # 9
             Line(p[20], p[21]),  # 10
             Line(p[22], p[23]),  # 11
             Line(p[1], p[9]),  # 12
             Line(p[3], p[11]),  # 13
             Line(p[5], p[13]),  # 14
             Line(p[7], p[15]),  # 15
             Line(p[9], p[17]),  # 16
             Line(p[11], p[19]),  # 17
             Line(p[13], p[21]),  # 18
             Line(p[15], p[23]),  # 19
             Line(p[0], p[8]),  # 20
             Line(p[6], p[14]),  # 21
             Line(p[8], p[16]),  # 22
             Line(p[10], p[18]),  # 23
             Line(p[12], p[20]),  # 24
             Line(p[14], p[22]),  # 25
             Line(p[3], p[1]),  # 26
             Line(p[5], p[3]),  # 27
             Line(p[5], p[7]),  # 28
             Line(p[11], p[9]),  # 29
             Line(p[13], p[11]),  # 30
             Line(p[13], p[15]),  # 31
             Line(p[19], p[17]),  # 32
             Line(p[21], p[19]),  # 33
             Line(p[21], p[23]),  # 34
             Line(p[2], p[0]),  # 35
             Line(p[4], p[6]),  # 36
             Line(p[10], p[8]),  # 37
             Line(p[12], p[14]),  # 38
             Line(p[4], p[2]),  # 39
             Line(p[12], p[10]),  # 40
             Line(p[18], p[16]),  # 41
             Line(p[20], p[18]),  # 42
             Line(p[20], p[22]),  # 43
             Line(p[24], p[30]),  # 44
             Line(p[25], p[31]),  # 45
             Line(p[26], p[32]),  # 46
             Line(p[27], p[33]),  # 47
             Line(p[28], p[34]),  # 48
             Line(p[29], p[35]),  # 49
             Line(p[26], p[28]),  # 50
             Line(p[27], p[29]),  # 51
             Line(p[30], p[32]),  # 52
             Line(p[31], p[33]),  # 53
             Line(p[32], p[34]),  # 54
             Line(p[33], p[35])  # 55
             ]

        c = [wsf.extract('U', 0), wsf.extract('U', 1)]
        c2, c3, c4 = Crv.split(outer_root, obrk_root)
        c5, c6, c7 = Crv.split(outer_tip, obrk_tip)
        c8, c9, c10 = Crv.split(outer_far, obrk_far)
        c11, c12, c13 = Crv.split(crv_root, brk_root)
        c14, c15, c16 = Crv.split(crv_tip, brk_tip)
        c17, c18, c19 = Crv.split(crv_far, brk_far)
        c20 = wsf.extract('U', brk_root[0])
        c21 = wsf.extract('U', brk_root[1])
        c22 = fsf.extract('U', brk_tip[0])
        c23 = fsf.extract('U', brk_tip[1])
        c3.reverse()
        c4.reverse()
        c6.reverse()
        c7.reverse()
        c9.reverse()
        c10.reverse()
        c12.reverse()
        c13.reverse()
        c15.reverse()
        c16.reverse()
        c18.reverse()
        c19.reverse()
        c.append(c2)
        c.append(c3)
        c.append(c4)
        c.append(c5)
        c.append(c6)
        c.append(c7)
        c.append(c8)
        c.append(c9)
        c.append(c10)
        c.append(c11)
        c.append(c12)
        c.append(c13)
        c.append(c14)
        c.append(c15)
        c.append(c16)
        c.append(c17)
        c.append(c18)
        c.append(c19)
        c.append(c20)
        c.append(c21)
        c.append(c22)
        c.append(c23)

        ts1 = Surf.split(wsf, brk_root, [])
        s0 = ts1[0][0]
        s1 = ts1[1][0]
        s2 = ts1[2][0]
        ts2 = Surf.split(fsf, brk_tip, [])
        s3 = ts2[0][0]
        s4 = ts2[1][0]
        s5 = ts2[2][0]
        s = [s0, s1, s2, s3, s4, s5]

        '''Node number and distribution'''
        n = np.full(8, enn, int)
        u0 = hyperbolic_tangent(n[0], 8)
        u1 = double_exponential(n[1], 0.5, 1.5, 0.5)
        u2 = uniform(n[2])
        u3 = single_exponential(n[3], 5)
        u4 = hyperbolic_tangent(n[4], 5)
        u5 = double_exponential(n[5], 0.5, 1.2, 0.5)
        u6 = double_exponential(n[6], 0.5, 1.5, 0.5)
        u7 = uniform(n[7])
        knot_dist = [u0, u1, u2, u3, u4, u5, u6, u7]

        '''Construct blocks'''
        blk_list = []
        blk_param_list = []

        b0_tfi_grid = LinearTFI3D.from_edges(l[1], l[26], l[0], l[35], l[5], l[29], l[4], l[37], c[0], l[13], l[12], l[20])
        blk_list.append(b0_tfi_grid)
        blk_param_list.append([knot_dist[3], knot_dist[0], knot_dist[2]])

        b1_tfi_grid = LinearTFI3D.from_edges(l[2], l[27], l[1], l[39], l[6], l[30], l[5], l[40], c[1], l[14], l[13], c[0])
        blk_list.append(b1_tfi_grid)
        blk_param_list.append([knot_dist[3], knot_dist[7], knot_dist[2]])

        b2_tfi_grid = LinearTFI3D.from_edges(l[36], l[3], l[28], l[2], l[38], l[7], l[31], l[6], c[1], l[21], l[15], l[14])
        blk_list.append(b2_tfi_grid)
        blk_param_list.append([knot_dist[0], knot_dist[3], knot_dist[2]])

        b3_tfi_grid = LinearTFI3D.from_edges(l[5], l[29], l[4], l[37], l[9], l[32], l[8], l[41], l[23], l[17], l[16], l[22])
        blk_list.append(b3_tfi_grid)
        blk_param_list.append([knot_dist[3], knot_dist[0], knot_dist[4]])

        b4_tfi_grid = LinearTFI3D.from_edges(l[6], l[30], l[5], l[40], l[10], l[33], l[9], l[42], l[24], l[18], l[17], l[23])
        blk_list.append(b4_tfi_grid)
        blk_param_list.append([knot_dist[3], knot_dist[7], knot_dist[4]])

        b5_tfi_grid = LinearTFI3D.from_edges(l[38], l[7], l[31], l[6], l[43], l[11], l[34], l[10], l[24], l[25], l[19], l[18])
        blk_list.append(b5_tfi_grid)
        blk_param_list.append([knot_dist[0], knot_dist[3], knot_dist[4]])

        b6_s1 = deepcopy(s[0])
        b6_s2 = LinearTFI2D(c[2], l[20], c[5], l[52])
        b6_s3 = LinearTFI2D(c[0], l[35], l[20], l[37])
        b6_s4 = LinearTFI2D(c[20], l[44], l[52], l[46])
        b6_s5 = LinearTFI2D(l[35], c[11], l[44], c[2])
        b6_s6 = LinearTFI2D(l[37], c[14], l[46], c[5])

        b6_tfi_grid = LinearTFI3D(lambda v, w: b6_s1(v, w),
                                  lambda v, w: b6_s2(v, w),
                                  lambda w, u: b6_s3(w, u),
                                  lambda w, u: b6_s4(w, u),
                                  lambda u, v: b6_s5(u, v),
                                  lambda u, v: b6_s6(u, v))

        blk_list.append(b6_tfi_grid)
        blk_param_list.append([knot_dist[0], knot_dist[1], knot_dist[2]])

        b7_s1 = LinearTFI2D(l[45], c[21], l[47], l[53])
        b7_s2 = LinearTFI2D(l[44], c[20], l[46], l[52])
        b7_s3 = deepcopy(s[1])
        b7_s3.reverse('U')
        b7_s3.swap()
        b7_s4 = LinearTFI2D(l[53], c[3], l[52], c[6])
        b7_s5 = LinearTFI2D(c[12], l[45], c[3], l[44])
        b7_s6 = LinearTFI2D(c[15], l[47], c[6], l[46])

        b7_tfi_grid = LinearTFI3D(lambda v, w: b7_s1(v, w),
                                  lambda v, w: b7_s2(v, w),
                                  lambda w, u: b7_s3(w, u),
                                  lambda w, u: b7_s4(w, u),
                                  lambda u, v: b7_s5(u, v),
                                  lambda u, v: b7_s6(u, v))

        blk_list.append(b7_tfi_grid)
        blk_param_list.append([knot_dist[5], knot_dist[0], knot_dist[2]])

        b8_s1 = LinearTFI2D(l[36], c[1], l[38], l[21])
        b8_s2 = LinearTFI2D(l[45], c[21], l[47], l[53])
        b8_s3 = deepcopy(s[2])
        b8_s3.reverse('U')
        b8_s3.swap()
        b8_s4 = LinearTFI2D(l[21], c[4], l[53], c[7])
        b8_s5 = LinearTFI2D(c[13], l[36], c[4], l[45])
        b8_s6 = LinearTFI2D(c[16], l[38], c[7], l[47])

        b8_tfi_grid = LinearTFI3D(lambda v, w: b8_s1(v, w),
                                  lambda v, w: b8_s2(v, w),
                                  lambda w, u: b8_s3(w, u),
                                  lambda w, u: b8_s4(w, u),
                                  lambda u, v: b8_s5(u, v),
                                  lambda u, v: b8_s6(u, v))

        blk_list.append(b8_tfi_grid)
        blk_param_list.append([knot_dist[6], knot_dist[0], knot_dist[2]])

        b9_s1 = deepcopy(s[3])
        b9_s2 = LinearTFI2D(c[5], l[22], c[8], l[54])
        b9_s3 = LinearTFI2D(l[23], l[37], l[22], l[41])
        b9_s4 = LinearTFI2D(c[22], l[46], l[54], l[48])
        b9_s5 = LinearTFI2D(l[37], c[14], l[46], c[5])
        b9_s6 = LinearTFI2D(l[41], c[17], l[48], c[8])

        b9_tfi_grid = LinearTFI3D(lambda v, w: b9_s1(v, w),
                                  lambda v, w: b9_s2(v, w),
                                  lambda w, u: b9_s3(w, u),
                                  lambda w, u: b9_s4(w, u),
                                  lambda u, v: b9_s5(u, v),
                                  lambda u, v: b9_s6(u, v))

        blk_list.append(b9_tfi_grid)
        blk_param_list.append([knot_dist[0], knot_dist[1], knot_dist[4]])

        b10_s1 = LinearTFI2D(l[47], c[23], l[49], l[55])
        b10_s2 = LinearTFI2D(l[46], c[22], l[48], l[54])
        b10_s3 = deepcopy(s[4])
        b10_s3.reverse('U')
        b10_s3.swap()
        b10_s4 = LinearTFI2D(l[55], c[6], l[54], c[9])
        b10_s5 = LinearTFI2D(c[15], l[47], c[6], l[46])
        b10_s6 = LinearTFI2D(c[18], l[49], c[9], l[48])

        b10_tfi_grid = LinearTFI3D(lambda v, w: b10_s1(v, w),
                                   lambda v, w: b10_s2(v, w),
                                   lambda w, u: b10_s3(w, u),
                                   lambda w, u: b10_s4(w, u),
                                   lambda u, v: b10_s5(u, v),
                                   lambda u, v: b10_s6(u, v))

        blk_list.append(b10_tfi_grid)
        blk_param_list.append([knot_dist[5], knot_dist[0], knot_dist[4]])

        b11_s1 = LinearTFI2D(l[38], l[24], l[43], l[25])
        b11_s2 = LinearTFI2D(l[47], c[23], l[49], l[55])
        b11_s3 = deepcopy(s[5])
        b11_s3.reverse('U')
        b11_s3.swap()
        b11_s4 = LinearTFI2D(l[25], c[7], l[55], c[10])
        b11_s5 = LinearTFI2D(c[16], l[38], c[7], l[47])
        b11_s6 = LinearTFI2D(c[19], l[43], c[10], l[49])

        b11_tfi_grid = LinearTFI3D(lambda v, w: b11_s1(v, w),
                                   lambda v, w: b11_s2(v, w),
                                   lambda w, u: b11_s3(w, u),
                                   lambda w, u: b11_s4(w, u),
                                   lambda u, v: b11_s5(u, v),
                                   lambda u, v: b11_s6(u, v))

        blk_list.append(b11_tfi_grid)
        blk_param_list.append([knot_dist[6], knot_dist[0], knot_dist[4]])

        b12_s1 = deepcopy(s[5])
        b12_s1.reverse('U')
        b12_s2 = deepcopy(s[3])
        b12_s3 = LinearTFI2D(l[24], l[40], l[23], l[42])
        b12_s4 = deepcopy(s[4])
        b12_s4.reverse('U')
        b12_s4.swap()
        b12_s5 = LinearTFI2D(l[40], c[16], c[15], c[14])
        b12_s6 = LinearTFI2D(l[42], c[19], c[18], c[17])

        b12_tfi_grid = LinearTFI3D(lambda v, w: b12_s1(v, w),
                                   lambda v, w: b12_s2(v, w),
                                   lambda w, u: b12_s3(w, u),
                                   lambda w, u: b12_s4(w, u),
                                   lambda u, v: b12_s5(u, v),
                                   lambda u, v: b12_s6(u, v))

        blk_list.append(b12_tfi_grid)
        blk_param_list.append([knot_dist[7], knot_dist[6], knot_dist[4]])

        def report(msg):
            print('Process {} : {}'.format(os.getpid(), msg))

        report('Calculating grid ...')
        for i in range(len(blk_list)):
            _u, _v, _w = blk_param_list[i]
            blk_list[i].calc_grid(_u, _v, _w)
            report('Block {} calculation done!'.format(i))

        '''网格, 边界条件, 邻接关系'''
        blk = [b0_tfi_grid.grid,
               b1_tfi_grid.grid,
               b2_tfi_grid.grid,
               b3_tfi_grid.grid,
               b4_tfi_grid.grid,
               b5_tfi_grid.grid,
               b6_tfi_grid.grid,
               b7_tfi_grid.grid,
               b8_tfi_grid.grid,
               b9_tfi_grid.grid,
               b10_tfi_grid.grid,
               b11_tfi_grid.grid,
               b12_tfi_grid.grid]

        bc = [(BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b0
              (BCType.Wall, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Symmetry, BCType.Interior),  # b1
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b2
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b3
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField),  # b4
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b5
              (BCType.Wall, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Symmetry, BCType.Interior),  # b6
              (BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b7
              (BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b8
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField),  # b9
              (BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b10
              (BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b11
              (BCType.Interior, BCType.Interior, BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField)  # b12
              ]

        adj = [((6, 3), (0, 1), 1, True),
               ((0, 2), (0, 0), 0, False),
               ((1, 4), (0, 3), 1, False),
               ((0, 4), (0, 0), 0, False),
               ((0, 0), (0, 5), 1, False),
               ((0, 6), (3, 5), 0, False),

               ((0, 0), (1, 1), 1, False),
               ((1, 2), (0, 0), 0, False),
               ((2, 1), (1, 3), 1, True),
               ((0, 0), (1, 5), 1, False),
               ((1, 6), (4, 5), 0, False),

               ((2, 2), (0, 0), 0, False),
               ((8, 1), (2, 3), 1, True),
               ((2, 4), (0, 0), 0, False),
               ((0, 0), (2, 5), 1, False),
               ((2, 6), (5, 5), 0, False),

               ((9, 3), (3, 1), 1, True),
               ((3, 2), (0, 0), 0, False),
               ((4, 4), (3, 3), 1, False),
               ((3, 4), (0, 0), 0, False),
               ((3, 6), (0, 0), 0, False),

               ((12, 3), (4, 1), 1, True),
               ((4, 2), (0, 0), 0, False),
               ((5, 1), (4, 3), 1, True),
               ((4, 6), (0, 0), 0, False),

               ((5, 2), (0, 0), 0, False),
               ((11, 1), (5, 3), 1, True),
               ((5, 4), (0, 0), 0, False),
               ((5, 6), (0, 0), 0, False),

               ((0, 0), (6, 1), 1, False),
               ((6, 2), (0, 0), 0, False),
               ((6, 4), (7, 2), 0, True),
               ((0, 0), (6, 5), 1, False),
               ((6, 6), (9, 5), 0, False),

               ((8, 2), (7, 1), 1, False),
               ((0, 0), (7, 3), 1, False),
               ((7, 4), (0, 0), 0, False),
               ((0, 0), (7, 5), 1, False),
               ((7, 6), (10, 5), 0, False),

               ((0, 0), (8, 3), 1, False),
               ((8, 4), (0, 0), 0, False),
               ((0, 0), (8, 5), 1, False),
               ((8, 6), (11, 5), 0, False),

               ((12, 2), (9, 1), 1, False),
               ((9, 2), (0, 0), 0, False),
               ((9, 4), (10, 2), 0, True),
               ((9, 6), (0, 0), 0, False),

               ((11, 2), (10, 1), 1, False),
               ((12, 4), (10, 3), 1, False),
               ((10, 4), (0, 0), 0, False),
               ((10, 6), (0, 0), 0, False),

               ((12, 1), (11, 3), 1, True),
               ((11, 4), (0, 0), 0, False),
               ((11, 6), (0, 0), 0, False),

               ((0, 0), (12, 5), 1, False),
               ((12, 6), (0, 0), 0, False)]

        '''构建MSH文件'''
        msh = XF_MSH.from_str3d_multi(blk, bc, adj)
        msh.save(fn)


class AirfoilTester(unittest.TestCase):
    def test_grid_gen(self):
        # airfoil, A, B, C, N0, N1, N2, N3
        data = [('SC(2)-0406', 30, 20, 50, 90, 60, 80, 3, 'none'),
                ('RAE2822', 30, 20, 50, 90, 60, 80, 3, 'none'),
                ('SC(2)-0406', 30, 20, 50, 90, 60, 80, 3, 'laplace'),
                ('RAE2822', 30, 20, 50, 90, 60, 80, 3, 'laplace'),
                ('NLF(1)-0414F', 30, 20, 50, 91, 61, 80, 3, 'thomas-middlecoff'),
                ('RAE2822', 30, 20, 50, 90, 60, 80, 3, 'thomas-middlecoff')]

        for k in range(len(data)):
            fn, la, lb, lc, n0, n1, n2, n3, smt = data[k]
            foil = Airfoil(fn)
            bunch = foil.gen_grid(la, lb, lc, n0, n1, n2, n3, leading_smooth=smt)
            p3d = bunch[1]
            p3d.save(fn + '_flowfield_grid-smooth={}.xyz'.format(smt))
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
