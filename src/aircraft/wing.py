import numpy as np
import os
from copy import deepcopy
import math
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy import interpolate
from src.iges.iges_core import IGES_Model
from src.nurbs.utility import equal, pnt_dist
from src.nurbs.curve import GlobalInterpolatedCrv
from src.nurbs.surface import Skinned
from settings import AIRFOIL_DIR


def update_airfoil_list():
    for f in os.listdir(AIRFOIL_DIR):
        base, ext = os.path.splitext(f)
        if ext == '.dat':
            Airfoil.AIRFOIL_LIST.append(base)


class Airfoil(object):
    AIRFOIL_LIST = []

    def __init__(self):
        """
        2D Airfoil, with chord length equals to 1.
        """

        self.name = None
        self.pts = None

    def __repr__(self):
        return self.name

    @property
    def front(self):
        total = self.pnt_num
        cx = self.pts[0][0]
        k = 1
        while k < total and self.pts[k][0] < cx:
            cx = self.pts[k][0]

        return self.pts[k - 1]

    @property
    def tail(self):
        return (self.pts[0] + self.pts[-1]) / 2

    @property
    def chord_len(self):
        return pnt_dist(self.front, self.tail)

    @property
    def is_blunt(self):
        return equal(norm(self.pts[0] - self.pts[-1], np.inf), 0.0)

    @property
    def pnt_num(self):
        return len(self.pts)

    def nurbs_rep(self, p=5, method='centripetal'):
        """
        NURBS Representation
        :param p: 插值次数
        :type p: int
        :param method: 插值参数化方法
        :type method: str
        :return: 翼型的NURBS全局插值曲线
        :rtype: GlobalInterpolatedCrv
        """

        return GlobalInterpolatedCrv(self.pts, p, method)

    @classmethod
    def read_pts(cls, fn):
        pts = []
        fin = open(os.path.join(AIRFOIL_DIR, fn + '.dat'))
        for line in fin:
            (x, y, z) = line.split()
            pts.append([x, y, z])
        fin.close()

        n = len(pts)
        ret = np.empty((n, 3), float)
        for k, pnt in enumerate(pts):
            for d in range(3):
                ret[k][d] = pnt[d]

        return ret

    @classmethod
    def from_file(cls, fn):
        """
        Construct airfoil from file.
        :param fn: File name
        :type fn: str
        :return: Airfoil Object
        :rtype: Airfoil
        """

        af = cls()
        af.name = fn
        af.pts = cls.read_pts(fn)
        return af

    def gen_msh(self):
        pass


class WingProfile(Airfoil):
    SEC_INTRINSIC_PARAM = ['Airfoil', 'Z(m)', 'X_front(m)', 'Y_front(m)', 'X_tail(m)', 'Y_tail(m)', 'Thickness Ratio']
    SEC_GEOM_PARAM = ['Airfoil', 'Z(m)', 'Length(m)', 'SweepBack(deg)', 'Twist(deg)', 'Dihedral(deg)', 'TwistPos', 'Thickness Ratio']

    def __init__(self, foil, ends, thickness_factor=1.0):
        """
        3D profile at certain position.
        :param foil: 翼型名称
        :type foil: str
        :param ends: 剖面起始端点
        :param thickness_factor: 翼型纵向拉伸系数
        :type thickness_factor: float
        """

        super(WingProfile, self).__init__()

        '''Inspect endings'''
        self.ending = np.copy(ends)
        if not equal(ends[0][2], ends[1][2]):
            raise AssertionError("Invalid ending coordinates in Z direction!")

        cl = self.chord_len
        if equal(cl, 0.0):
            raise ZeroDivisionError("Invalid ending coordinates in XY direction!")

        rotation = complex((ends[1][0] - ends[0][0]) / cl, (ends[1][1] - ends[0][1]) / cl)

        '''Build section'''
        self.name = foil
        self.pts = np.copy(super(WingProfile, self).read_pts(foil))
        n = self.pnt_num
        for i in range(n):
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
        ending = np.empty((2, 3), float)
        ending[0] = front
        ending[1] = tail

        return cls(foil, ending, thickness_factor)


class WingFrame(object):
    def __init__(self, xf, xt, yf, yt, z):
        """
        机翼外框描述
        :param xf: 前缘x坐标的参数方程
        :param xt: 后缘x坐标的参数方程
        :param yf: 前缘y坐标的参数方程
        :param yt: 后缘y坐标的参数方程
        :param z: z坐标的参数方程
        """

        self.f_xf = deepcopy(xf)
        self.f_xt = deepcopy(xt)
        self.f_yf = deepcopy(yf)
        self.f_yt = deepcopy(yt)
        self.f_z = deepcopy(z)

    def x_front(self, u):
        return self.f_xf(u)

    def x_tail(self, u):
        return self.f_xt(u)

    def y_front(self, u):
        return self.f_yf(u)

    def y_tail(self, u):
        return self.f_yt(u)

    def z(self, u):
        return self.f_z(u)

    def show(self, n=1000):
        u_dist = np.linspace(0, 1.0, n)
        z = np.empty(n, float)
        xf = np.empty(n, float)
        xt = np.empty(n, float)
        for k in range(n):
            z[k] = self.z(u_dist[k])
            xf[k] = self.x_front(u_dist[k])
            xt[k] = self.x_tail(u_dist[k])

        plt.figure()
        plt.plot(z, xf, label='Front')
        plt.plot(z, xt, label='Tail')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.show()


class BWBFrame(WingFrame):
    def __init__(self, c_root, c_mid, c_tip, b_mid, b_tip, alpha_root, alpha_mid, alpha_tip):
        """
        BWB(Blended Wing Body)构型飞行器参数化描述
        :param c_root: 根部弦长
        :type c_root: float
        :param c_mid: 中间弦长
        :type c_mid: float
        :param c_tip: 翼尖弦长
        :type c_tip: float
        :param b_mid: 内翼宽度
        :type b_mid: float
        :param b_tip: 机翼半展长
        :type b_tip: float
        :param alpha_root: 内翼平均后掠角
        :type alpha_root: float
        :param alpha_mid: 中段平均后掠角
        :type alpha_mid: float
        :param alpha_tip: 翼尖平均后掠角
        :type alpha_tip: float
        """

        self.Cr = c_root
        self.Cm = c_mid
        self.Ct = c_tip
        self.Bm = b_mid
        self.Bt = b_tip
        self.Ar = alpha_root
        self.Am = alpha_mid
        self.At = alpha_tip

        '''Calculate pivots on each curve'''
        front_pnt = np.empty((3, 3), float)
        tail_pnt = np.empty((3, 3), float)
        front_pnt[0] = np.zeros(3)
        tail_pnt[0] = np.array([self.Cr, 0, 0])
        front_pnt[1] = np.array([self.Bm * math.tan(math.radians(self.Am)), 0, self.Bm])
        tail_pnt[1] = np.array([front_pnt[1][0] + self.Cm, 0, self.Bm])
        front_pnt[2] = np.array([front_pnt[1][0] + (self.Bt - self.Bm) * math.tan(math.radians(self.At)), 0, self.Bt])
        tail_pnt[2] = np.array([front_pnt[2][0] + self.Ct, 0, self.Bt])

        '''Build interpolated functions'''
        u = np.array([0, self.Bm / self.Bt, 1.0])
        z = np.array([0, self.Bm, self.Bt])
        xf = np.array([front_pnt[0][0], front_pnt[1][0], front_pnt[2][0]])
        yf = np.array([front_pnt[0][1], front_pnt[1][1], front_pnt[2][1]])
        xt = np.array([tail_pnt[0][0], tail_pnt[1][0], tail_pnt[2][0]])
        yt = np.array([tail_pnt[0][1], tail_pnt[1][1], tail_pnt[2][1]])

        super(BWBFrame, self).__init__(interpolate.make_interp_spline(u, xf, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       interpolate.make_interp_spline(u, xt, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       interpolate.make_interp_spline(u, yf, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       interpolate.make_interp_spline(u, yt, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       lambda t: z[1] * t / u[1] if t <= u[1] else z[1] + (z[2] - z[1]) * (t - u[1]) / (u[2] - u[1]))


class Wing(object):
    def __init__(self, section_list):
        """
        从剖面序列构造机翼
        :param section_list: 机翼剖面序列
        """

        self.section = []
        for elem in section_list:
            self.section.append(deepcopy(elem))

    @property
    def section_num(self):
        return len(self.section)

    @property
    def root(self):
        return self.section[0].nurbs_rep()

    @property
    def tip(self):
        return self.section[-1].nurbs_rep()

    def front(self, q=3, method='chord'):
        n = self.section_num
        front_pts = np.zeros((n, 3))
        for i in range(n):
            front_pts[i] = self.section[i].front

        return GlobalInterpolatedCrv(front_pts, q, method)

    def surf(self, p=3, q=3):
        """
        构建机翼轮廓曲线、曲面
        :param p: U方向次数
        :type p: int
        :param q: V方向次数
        :type q: int
        :return: 机翼蒙皮曲面
        :rtype: Skinned
        """

        profile_list = []
        for elem in self.section:
            profile_list.append(elem.nurbs_rep(p))
        return Skinned(profile_list, p, q)

    @property
    def tail_up(self):
        sk = self.surf()
        return sk.extract('U', 0)

    @property
    def tail_down(self):
        sk = self.surf()
        return sk.extract('U', 1)

    def iges_model(self, fn, p=3, q=3, mirror=True):
        """
        生成机翼相应的IGES模型
        :param fn: 文件名
        :type fn: str
        :param p: U方向次数
        :type p: int
        :param q: V方向次数
        :type q: int
        :param mirror: 是否生成对称部分
        :type mirror: bool
        :return: 可用于后续生成IGS文件的IGES_Model对象
        :rtype: IGES_Model
        """

        wing_model = IGES_Model(fn)

        '''前缘曲线'''
        wing_model.add_entity(self.front(q).to_iges())

        '''剖面'''
        for elem in self.section:
            wing_model.add_entity(elem.nurbs_rep(p).to_iges())

        '''蒙皮'''
        sk = self.surf(p, q)
        wing_model.add_entity(sk.to_iges())

        '''镜像'''
        if mirror:
            msk = deepcopy(sk)
            for i in range(msk.n + 1):
                for j in range(msk.m + 1):
                    msk.Pw[i][j][2] *= -1

            wing_model.add_entity(msk.to_iges())

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
    def from_frame(cls):
        pass
