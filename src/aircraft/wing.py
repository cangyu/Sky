import math
import os
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
from numpy.linalg import norm
from src.iges.iges_core import IGES_Model
from src.nurbs.curve import GlobalInterpolatedCrv
from src.nurbs.surface import Skinned
from src.nurbs.utility import equal, pnt_dist
from src.aircraft.frame import WingFrame
from settings import AIRFOIL_DIR


class Airfoil(object):
    AIRFOIL_LIST = []

    def __init__(self):
        """
        2D Airfoil, with chord length equals to 1.
        """

        self.name = None
        self.pts = None

    def __repr__(self):
        ret = "Airfoil: {}\nNbrOfPnt: {}\nCoordinates(from tail-up to tail-down):\n".format(self.name, len(self.pts))
        crd = []
        for k, pnt in enumerate(self.pts):
            crd.append("{}: ({}, {}, {})".format(k + 1, pnt[0], pnt[1], pnt[2]))

        ret += "\n".join(crd)
        return ret

    def __str__(self):
        return self.name

    @classmethod
    def update_airfoil_list(cls):
        for f in os.listdir(AIRFOIL_DIR):
            base, ext = os.path.splitext(f)
            if ext == '.dat':
                cls.AIRFOIL_LIST.append(base)

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
    def tail_up(self):
        return self.pts[0]

    @property
    def tail_down(self):
        return self.pts[-1]

    @property
    def chord_len(self):
        return pnt_dist(self.front, self.tail)

    @property
    def is_blunt(self):
        return equal(norm(self.pts[0] - self.pts[-1], np.inf), 0.0)

    @property
    def pnt_num(self):
        return len(self.pts)

    def show(self):
        (px, py, pz) = zip(*self.pts)
        plt.plot(px, py)
        plt.gca().set_aspect('equal')
        plt.show()

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
    def from_database(cls, fn):
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


if __name__ == '__main__':
    af = Airfoil.from_database('M6')
    print(af)
    print(repr(af))
