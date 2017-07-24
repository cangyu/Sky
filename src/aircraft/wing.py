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
from src.nurbs.surface import GlobalInterpolatedSurf, Skinned
from settings import AIRFOIL_DIR

BWB_SEC_PARAM = ['Airfoil', 'Thickness Ratio', 'Z(m)', 'X_front(m)', 'Y_front(m)', 'X_tail(m)', 'Y_tail(m)']
AIRFOIL_LIST = []


def update_airfoil_list():
    for f in os.listdir(AIRFOIL_DIR):
        base, ext = os.path.splitext(f)
        if ext == '.dat':
            AIRFOIL_LIST.append(base)


class Airfoil(object):
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

    @property
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

        return pts

    @classmethod
    def from_file(cls, fn):
        """
        Construct airfoil from file.
        :param fn: File name
        :type fn: str
        :return: Airfoil Object
        :rtype: Airfoil
        """

        pts = cls.read_pts(fn)
        af = cls()
        af.name = fn
        af.pts = np.copy(pts)
        return af


class WingProfile(Airfoil):
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
    def __init__(self, airfoil_list, thickness_list, z_list, xf_list, yf_list, xt_list, yt_list):
        self.airfoil = airfoil_list
        self.thickness = thickness_list
        self.z = z_list
        self.xf = xf_list
        self.yf = yf_list
        self.xt = xt_list
        self.yt = yt_list

        self.surf = None
        self.root = None
        self.tip = None
        self.front = None
        self.tail_up = None
        self.tail_down = None

    @property
    def section_num(self):
        return len(self.airfoil)

    def build_sketch(self, p=5, q=5):
        """
        构建机翼轮廓曲线、曲面
        :param p: U方向次数
        :param q: V方向次数
        :return: None
        """

        '''剖面'''
        profile = []
        for i in range(0, self.n):
            epts = np.array([[self.xf[i], self.yf[i], self.z[i]],
                             [self.xt[i], self.yt[i], self.z[i]]])
            wp = WingProfile(self.airfoil[i], epts, self.thickness[i])
            profile.append(wp)

        self.root = profile[0].nurbs_rep
        self.tip = profile[-1].nurbs_rep

        '''曲面全局插值'''
        nn = len(profile[0].pts)
        mm = self.n
        surf_pts = np.zeros((nn, mm, 3))
        for i in range(0, nn):
            for j in range(0, mm):
                surf_pts[i][j] = np.copy(profile[j].pts[i])
        self.surf = GlobalInterpolatedSurf(surf_pts, p, q)

        '''前后缘采样点'''
        front_pts = np.zeros((self.n, 3))
        tail_up_pts = np.zeros((self.n, 3))
        tail_down_pts = np.zeros((self.n, 3))
        for i in range(self.n):
            front_pts[i][0] = self.xf[i]
            front_pts[i][1] = self.yf[i]
            front_pts[i][2] = self.z[i]
            tail_up_pts[i] = np.copy(profile[i].pts[0])
            tail_down_pts[i] = np.copy(profile[i].pts[-1])

        self.front = GlobalInterpolatedCrv(front_pts, q, 'chord')
        self.tail_up = self.surf.extract('U', 0)
        self.tail_down = self.surf.extract('U', 1)

    @property
    def geom(self):
        return self.surf, self.root, self.tip, self.front, self.tail_up, self.tail_down

    def write(self, fn, p=5, q=5, mirror=True):
        wing_model = IGES_Model(fn)

        '''前后缘采样点'''
        front_pts = np.zeros((self.n, 3))
        tail_pts = np.zeros((self.n, 3))
        for i in range(0, self.n):
            front_pts[i][0] = self.xf[i]
            front_pts[i][1] = self.yf[i]
            tail_pts[i][0] = self.xt[i]
            tail_pts[i][1] = self.yt[i]
            tail_pts[i][2] = front_pts[i][2] = self.z[i]

        '''
        for i in range(0, self.n):
            wing_model.add_entity(IGES_Entity116(self.xf[i], self.yf[i], self.z[i]))
            wing_model.add_entity(IGES_Entity116(self.xt[i], self.yt[i], self.z[i]))
        '''

        '''前后缘曲线'''
        front_crv = GlobalInterpolatedCrv(front_pts, 5, 'chord')
        tail_crv = GlobalInterpolatedCrv(tail_pts, 5, 'chord')
        wing_model.add_entity(front_crv.to_iges(0, 0, [0, 0, 0]))
        wing_model.add_entity(tail_crv.to_iges(0, 0, [0, 0, 0]))

        '''根梢弦线'''
        # wing_model.add_entity(IGES_Entity110([self.xf[0], self.yf[0], self.z[0]], [self.xt[0], self.yt[0], self.z[0]]))
        # wing_model.add_entity(IGES_Entity110([self.xf[- 1], self.yf[- 1], self.z[- 1]], [self.xt[- 1], self.yt[- 1], self.z[- 1]]))

        '''剖面'''
        profile = []
        for i in range(0, self.n):
            epts = np.array([[self.xf[i], self.yf[i], self.z[i]],
                             [self.xt[i], self.yt[i], self.z[i]]])
            wp = WingProfile(self.airfoil[i], epts, self.thickness[i]).nurbs_rep
            profile.append(wp)
            wing_model.add_entity(wp.to_iges(1, 0, [0, 0, 1]))
            # add_pnt(wing_model, wp.pts[0])
            # add_pnt(wing_model, wp.pts[-1])

        '''全局插值'''
        self.surf = Skinned(profile, p, q)
        wing_model.add_entity(self.surf.to_iges())

        if mirror:
            msrf = deepcopy(self.surf)
            for i in range(msrf.n + 1):
                for j in range(msrf.m + 1):
                    msrf.Pw[i][j][2] *= -1

            wing_model.add_entity(msrf.to_iges())

        return wing_model
