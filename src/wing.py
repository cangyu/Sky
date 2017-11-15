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
from misc import pnt_dist, read_airfoil_pts, pnt_pan
from spacing import hyperbolic_tangent, uniform, single_exponential, double_exponential
from iges import Model, Entity116, Entity110
from nurbs import Crv, ConicArc, Line, Spline, ConicArc, Surf, Skinned, RuledSurf
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
        :param kwargs: Extra options on smoothing and spacing.
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

    def gen_grid(self, *args, **kwargs):
        """
        Generate the multi-block grid for a simplified wing.
        :param args: Containing the geometric description and node distribution of the grid.
        :param kwargs: Extra options on smoothing and spacing.
        :return: The wire-frame, plot3d-grid and fluent-grid(with predefined BC) of the flow field.
        """

        assert len(args) == 5

        a, b, c, d, n = args
        assert len(n) == 8

        '''Node distribution'''
        u = [hyperbolic_tangent(n[0], 8),  # u0
             double_exponential(n[1], 0.5, 1.5, 0.5),  # u1
             uniform(n[2]),  # u2
             single_exponential(n[3], 5),  # u3
             hyperbolic_tangent(n[4], 5),  # u4
             double_exponential(n[5], 0.5, 1.2, 0.5),  # u5
             double_exponential(n[6], 0.5, 1.5, 0.5),  # u6
             uniform(n[7])]  # u7

        '''
        wsf: Wing surface.
        rt: Wing root curve.
        rl: Wing root chord length.
        tp: Wing tip curve.
        z: Z-Dim position of root, tip and far.
        spn: Wing span(From root to tip).
        far: Wing tip extension on far-field.
        wing_brk: Splitting position on wing surf in U direction to construct separated blocks.
        inlet_brk: Splitting position on inlet surf in U direction to construct separated blocks.
        '''
        wsf = self.surf
        rt = wsf.extract('V', 0)
        rl = self.profile[0].chord_len
        tp = wsf.extract('V', 1)
        z = np.array([self.profile[0].z, self.profile[-1].z, self.profile[-1].z + d])
        spn = z[1] - z[0]
        far = WingProfile.from_geom_param(self.profile[-1].name, z[2], rl, 0, 0, 0, thickness_factor=2).crv
        fsf = RuledSurf(tp, far)
        wing_brk = kwargs['wing_brk_param'] if 'wing_brk_param' in kwargs else (0.45, 0.55)
        inlet_brk = kwargs['inlet_brk'] if 'inlet_brk' in kwargs else (0.3, 0.72)

        assert len(wing_brk) == 2
        assert len(inlet_brk) == 2

        t1 = (rl, b, z[0])
        t2 = (rl, -b, z[0])
        inlet1 = ConicArc(t1, (-1, 0, 0), t2, (1, 0, 0), (rl - a, 0, z[0]))
        t1 = pnt_pan(inlet1.start, (0, 0, spn))
        t2 = pnt_pan(inlet1.end, (0, 0, spn))
        inlet2 = ConicArc(t1, (-1, 0, 0), t2, (1, 0, 0), (rl - a, 0, z[1]))
        t1 = pnt_pan(inlet2.start, (0, 0, d))
        t2 = pnt_pan(inlet2.end, (0, 0, d))
        inlet3 = ConicArc(t1, (-1, 0, 0), t2, (1, 0, 0), (rl - a, 0, z[2]))

        '''Points, lines, curves, surfs'''
        p = np.zeros((36, 3))
        p[0] = inlet1.start
        p[1] = (p[0][0] + c, p[0][1], z[0])
        p[2] = wsf(0, 0)
        p[3] = (p[1][0], p[2][1], z[0])
        p[4] = wsf(1, 0)
        p[5] = (p[3][0], p[4][1], z[0])
        p[6] = inlet1.end
        p[7] = (p[5][0], p[6][1], z[0])
        p[8] = inlet2.start
        p[9] = pnt_pan(p[1], (0, 0, spn))
        p[10] = wsf(0, 1)
        p[11] = (p[9][0], p[10][1], z[1])
        p[12] = wsf(1, 1)
        p[13] = (p[11][0], p[12][1], z[1])
        p[14] = inlet2.end
        p[15] = pnt_pan(p[7], (0, 0, spn))
        p[16] = inlet3.start
        p[17] = pnt_pan(p[9], (0, 0, d))
        p[18] = far.start
        p[19] = (p[17][0], p[18][1], z[2])
        p[20] = far.end
        p[21] = (p[19][0], p[20][1], z[2])
        p[22] = inlet3.end
        p[23] = pnt_pan(p[15], (0, 0, d))
        p[24] = rt(wing_brk[0])
        p[25] = rt(wing_brk[1])
        p[26] = tp(wing_brk[0])
        p[27] = tp(wing_brk[1])
        p[28] = far(wing_brk[0])
        p[29] = far(wing_brk[1])
        p[30] = inlet1(inlet_brk[0])
        p[31] = inlet1(inlet_brk[1])
        p[32] = inlet2(inlet_brk[0])
        p[33] = inlet2(inlet_brk[1])
        p[34] = inlet3(inlet_brk[0])
        p[35] = inlet3(inlet_brk[1])

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
             Line(p[33], p[35])]  # 55

        c = [wsf.extract('U', 0), wsf.extract('U', 1)]  # c0, c1
        c += Crv.split(inlet1, inlet_brk)  # c2, c3, c4
        c += Crv.split(inlet2, inlet_brk)  # c5, c6, c7
        c += Crv.split(inlet3, inlet_brk)  # c8, c9, c10
        c += Crv.split(rt, wing_brk)  # c11, c12, c13
        c += Crv.split(tp, wing_brk)  # c14, c15, c16
        c += Crv.split(far, wing_brk)  # c17, c18, c19
        c += [wsf.extract('U', bp) for bp in wing_brk]  # c20, c21
        c += [fsf.extract('U', bp) for bp in wing_brk]  # c22, c23
        for k in (3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19):
            c[k].reverse()

        ts1 = Surf.split(wsf, wing_brk, [])
        ts2 = Surf.split(fsf, wing_brk, [])
        s = [ts1[0][0], ts1[1][0], ts1[2][0], ts2[0][0], ts2[1][0], ts2[2][0]]

        '''IGES Model'''
        iges_model = Model()
        for pnt in p:
            iges_model.add(Entity116(pnt[0], pnt[1], pnt[2]))
        for line in l:
            iges_model.add(Entity110(line.start, line.end))
        for crv in c:
            iges_model.add(crv.to_iges())
        for surf in s:
            iges_model.add(surf.to_iges())

        return iges_model


class AirfoilTestCase(unittest.TestCase):
    def test_2d_grid(self):
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

    def test_3d_grid(self):
        foil = ['SC(2)-0414', 'SC(2)-0414', 'SC(2)-0612', 'SC(2)-0712', 'SC(2)-0710', 'SC(2)-0710', 'SC(2)-0710', 'SC(2)-1010', 'SC(2)-1010', 'SC(2)-1006', 'SC(2)-0706', 'SC(2)-0706', 'SC(2)-0606', 'SC(2)-0406']
        length = 1.0 + uniform(len(foil))[::-1] * 3
        thickness_factor = np.ones(len(foil))
        z_offset = uniform(len(foil)) * 20
        sweep_back = np.full(len(foil), 25.0)
        twist = np.zeros(len(foil))
        dihedral = np.zeros(len(foil))
        twist_pos = np.full(len(foil), 1)
        y_ref = np.zeros(len(foil))
        wing = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)
        iges_model = wing.gen_grid(80, 20, 120, 300, (60, 40, 200, 100, 90, 80, 40, 10))
        iges_model.save('test.igs')
        self.assertTrue(True)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(AirfoilTestCase('test_3d_grid'))
    runner = unittest.TextTestRunner()
    runner.run(suite)
