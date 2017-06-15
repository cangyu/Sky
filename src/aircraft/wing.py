import numpy as np
import os
from src.iges.iges_core import IGES_Model
from src.iges.iges_entity110 import IGES_Entity110
from src.iges.iges_entity116 import IGES_Entity116
from src.nurbs.utility import equal, pnt_dist
from src.nurbs.curve import GlobalInterpolatedCrv
from src.nurbs.surface import GlobalInterpolatedSurf

BWB_SEC_PARAM = ['Airfoil', 'Thickness Ratio', 'Z(m)', 'X_front(m)', 'Y_front(m)', 'X_tail(m)', 'Y_tail(m)']

AIRFOIL_DIR = '../airfoil/'
AIRFOIL_LIST = []


def update_airfoil_list():
    for f in os.listdir(AIRFOIL_DIR):
        base, ext = os.path.splitext(f)
        if ext == '.dat':
            AIRFOIL_LIST.append(base)


class Airfoil(object):
    def __init__(self, name):
        """
        2D Airfoil, with chord length equals to 1.
        """

        '''Read input data'''
        pts = []
        fin = open("{}{}.dat".format(AIRFOIL_DIR, name))
        for line in fin:
            (_x, _y, _z) = line.split()
            pts.append(np.array([_x, _y, _z]))
        fin.close()

        '''Reconstruct'''
        self.name = name
        self.n = len(pts)
        self.pts = np.zeros((self.n, 3))
        for i in range(0, self.n):
            self.pts[i] = np.copy(pts[i])


class WingProfile(object):
    def __init__(self, airfoil, ends, thickness_factor=1.0, p=5):
        """
        3D profile at given position.
        """

        self.airfoil = Airfoil(airfoil)
        self.ends = ends
        self.thickness = thickness_factor
        self.n = self.airfoil.n
        self.pts = np.copy(self.airfoil.pts)

        if not equal(ends[0][2], ends[1][2]):
            raise ValueError("Invalid ending coordinates in Z direction!")

        chordLen = pnt_dist(ends[0], ends[1])
        if equal(chordLen, 0.0):
            raise ValueError("Invalid ending coordinates in XY direction!")

        rotation = complex((ends[1][0] - ends[0][0]) / chordLen, (ends[1][1] - ends[0][1]) / chordLen)

        '''Stretch, Z offset and Thickness'''
        for i in range(0, self.n):
            self.pts[i][0] *= chordLen
            self.pts[i][1] *= (chordLen * thickness_factor)
            self.pts[i][2] = ends[0][2]

        '''Rotate around ends[0]'''
        for i in range(0, self.n):
            origin_vector = complex(self.pts[i][0], self.pts[i][1])
            origin_vector *= rotation

            self.pts[i][0] = origin_vector.real
            self.pts[i][1] = origin_vector.imag

        '''Move to ends[0]'''
        for i in range(0, self.n):
            self.pts[i][0] += ends[0][0]
            self.pts[i][1] += ends[0][1]

        '''NURBS Representation'''
        self.nurbs_rep = GlobalInterpolatedCrv(self.pts, p, 'centripetal')


def add_line(container: IGES_Model, pts, i: int, j: int):
    container.add_entity(IGES_Entity110([pts[i][0], pts[i][1], pts[i][2]],
                                        [pts[j][0], pts[j][1], pts[j][2]]))


def add_pnt(container: IGES_Model, pnt):
    container.add_entity(IGES_Entity116(pnt[0], pnt[1], pnt[2]))


class Wing(object):
    def __init__(self, airfoil_list, thickness_list, z_list, xf_list, yf_list, xt_list, yt_list):
        self.n = len(airfoil_list)
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
        self.tail_up = GlobalInterpolatedCrv(tail_up_pts, q, 'chord')
        self.tail_down = GlobalInterpolatedCrv(tail_down_pts, q, 'chord')

    @property
    def geom(self):
        return self.surf, self.root, self.tip, self.front, self.tail_up, self.tail_down

    def write(self, fn, p=5, q=5, mirror=True, farfield_box=False, H=80, L=600, W=250):
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

        for i in range(0, self.n):
            wing_model.add_entity(IGES_Entity116(self.xf[i], self.yf[i], self.z[i]))
            wing_model.add_entity(IGES_Entity116(self.xt[i], self.yt[i], self.z[i]))

        '''前后缘曲线'''
        front_crv = GlobalInterpolatedCrv(front_pts, 5, 'chord')
        tail_crv = GlobalInterpolatedCrv(tail_pts, 5, 'chord')
        wing_model.add_entity(front_crv.to_iges(0, 0, [0, 0, 0]))
        wing_model.add_entity(tail_crv.to_iges(0, 0, [0, 0, 0]))

        '''根梢弦线'''
        wing_model.add_entity(IGES_Entity110([self.xf[0], self.yf[0], self.z[0]], [self.xt[0], self.yt[0], self.z[0]]))
        wing_model.add_entity(IGES_Entity110([self.xf[- 1], self.yf[- 1], self.z[- 1]], [self.xt[- 1], self.yt[- 1], self.z[- 1]]))

        '''剖面'''
        profile = []
        for i in range(0, self.n):
            epts = np.array([[self.xf[i], self.yf[i], self.z[i]],
                             [self.xt[i], self.yt[i], self.z[i]]])
            wp = WingProfile(self.airfoil[i], epts, self.thickness[i])
            wing_model.add_entity(wp.nurbs_rep.to_iges(1, 0, [0, 0, 1]))
            add_pnt(wing_model, wp.pts[0])
            add_pnt(wing_model, wp.pts[-1])
            profile.append(wp)

        '''全局插值'''
        nn = len(profile[0].pts)
        mm = self.n
        surf_pts = np.zeros((nn, mm, 3))
        for i in range(0, nn):
            for j in range(0, mm):
                surf_pts[i][j] = np.copy(profile[j].pts[i])
        self.surf = GlobalInterpolatedSurf(surf_pts, p, q)
        wing_model.add_entity(self.surf.to_iges(0, 0, 0, 0))

        if mirror:
            for i in range(nn):
                for j in range(mm):
                    surf_pts[i][j][2] = -surf_pts[i][j][2]
            surf = GlobalInterpolatedSurf(surf_pts, p, q)
            wing_model.add_entity(surf.to_iges(0, 0, 0, 0))

        '''远场边框'''
        if farfield_box:
            farfield_pts = np.array([[-L / 2, -H / 2, 0],
                                     [-L / 2, H / 2, 0],
                                     [L / 2, H / 2, 0],
                                     [L / 2, -H / 2, 0],
                                     [-L / 2, -H / 2, W],
                                     [-L / 2, H / 2, W],
                                     [L / 2, H / 2, W],
                                     [L / 2, -H / 2, W]])

            for k in range(0, 8):
                wing_model.add_entity(IGES_Entity116(farfield_pts[k][0], farfield_pts[k][1], farfield_pts[k][2]))

            for k in range(0, 4):
                add_line(wing_model, farfield_pts, k, (k + 1) % 4)
                add_line(wing_model, farfield_pts, k + 4, (k + 5) % 4 + 4)
                add_line(wing_model, farfield_pts, k, k + 4)

        wing_model.write()
