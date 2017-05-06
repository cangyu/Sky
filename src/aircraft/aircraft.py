from src.iges.iges_entity110 import *
from src.iges.iges_entity116 import *
from src.iges.iges_entity126 import *
from src.aircraft.wing import *
from src.nurbs.curve import *
from src.nurbs.surface import *


class Aircraft(object):
    def __init__(self, airfoil_list, thickness_list, z_list, xf_list, yf_list, xt_list, yt_list):
        self.n = len(airfoil_list)
        self.airfoil = airfoil_list
        self.thickness = thickness_list
        self.z = z_list
        self.xf = xf_list
        self.yf = yf_list
        self.xt = xt_list
        self.yt = yt_list

    def generate(self, show_airfoil=True, show_surf=False, show_airfoil_cross_lines=False, show_airfoil_pts=False):
        plane = IGES_Model()

        # 前后缘采样点
        front_pts = np.zeros((self.n, 3), float)
        tail_pts = np.zeros((self.n, 3), float)
        for i in range(0, self.n):
            front_pts[i][0] = self.xf[i]
            front_pts[i][1] = self.yf[i]
            tail_pts[i][0] = self.xt[i]
            tail_pts[i][1] = self.yt[i]
            tail_pts[i][2] = front_pts[i][2] = self.z[i]
        front_crv = Curve(front_pts)
        tail_crv = Curve(tail_pts)

        for i in range(0, self.n):
            plane.AddPart(IGES_Entity116(self.xf[i], self.yf[i], self.z[i]))
            plane.AddPart(IGES_Entity116(self.xt[i], self.yt[i], self.z[i]))

        # 前缘曲线
        a, b, c = front_crv.generate()
        plane.AddPart(IGES_Entity126(front_crv.p, front_crv.n, 0, 0, 1, 0, a, b, c, 0.0, 1.0, np.array([0, 0, 0], float)))

        # 尾缘曲线
        a, b, c = tail_crv.generate()
        plane.AddPart(IGES_Entity126(tail_crv.p, tail_crv.n, 0, 0, 1, 0, a, b, c, 0.0, 1.0, np.array([0, 0, 0], float)))

        # 根梢弦线
        plane.AddPart(IGES_Entity110([self.xf[0], self.yf[0], self.z[0]],
                                     [self.xt[0], self.yt[0], self.z[0]]))
        plane.AddPart(IGES_Entity110([self.xf[self.n - 1], self.yf[self.n - 1], self.z[self.n - 1]],
                                     [self.xt[self.n - 1], self.yt[self.n - 1], self.z[self.n - 1]]))

        profile_pts = []
        profile_crv = []

        for i in range(0, self.n):
            epts = np.array([[self.xf[i], self.yf[i], self.z[i]], [self.xt[i], self.yt[i], self.z[i]]])
            wp = Wing_Profile(self.airfoil[i], epts, self.thickness[i])

            # 翼型离散点
            pts = wp.getPointList()
            profile_pts.append(pts)
            if show_airfoil_pts:
                for k in range(0, len(pts)):
                    plane.AddPart(IGES_Entity116(pts[k][0], pts[k][1], pts[k][2]))

            # 翼型NURBS曲线，默认5阶
            cc = Curve(pts)
            a, b, c = cc.generate()
            if show_airfoil:
                plane.AddPart(IGES_Entity126(cc.p, cc.n, 1, 0, 1, 0, a, b, c, 0.0, 1.0, np.array([0, 0, 1.0], float)))
                plane.AddPart(IGES_Entity110([pts[0][0], pts[0][1], pts[0][2]],
                                             [pts[-1][0], pts[-1][1], pts[-1][2]]))
            profile_crv.append(cc)

        if show_surf:
            surf = Skining(profile_crv)
            plane.AddPart(surf.generate())

        if show_airfoil_cross_lines:
            n = len(profile_pts[0])
            m = self.n
            for i in range(0, n):
                pts = np.zeros((m, 3), float)
                for j in range(0, m):
                    pts[j] = profile_pts[j][i]

                cc = Curve(pts)
                a, b, c = cc.generate('chord')
                plane.AddPart(IGES_Entity126(cc.p, cc.n, 0, 0, 1, 0, a, b, c, 0.0, 1.0, np.array([0, 0, 0], float)))

        # 尾缘
        tail_up_pts = []
        tail_down_pts = []
        for i in range(0, len(profile_pts)):
            tail_up_pts.append(profile_pts[i][0])
            tail_down_pts.append(profile_pts[i][-1])

        tail_up_crv = Curve(tail_up_pts)
        tail_down_crv = Curve(tail_down_pts)

        a, b, c = tail_up_crv.generate()
        plane.AddPart(
            IGES_Entity126(tail_up_crv.p, tail_up_crv.n, 0, 0, 1, 0, a, b, c, 0.0, 1.0, np.array([0, 0, 0], float)))

        a, b, c = tail_down_crv.generate()
        plane.AddPart(
            IGES_Entity126(tail_down_crv.p, tail_down_crv.n, 0, 0, 1, 0, a, b, c, 0.0, 1.0, np.array([0, 0, 0], float)))

        H = 80
        L = 600
        W = 250
        # 远场边框
        farfield_pts = np.array([[-L / 2, -H / 2, 0],
                                 [-L / 2, H / 2, 0],
                                 [L / 2, H / 2, 0],
                                 [L / 2, -H / 2, 0],
                                 [-L / 2, -H / 2, W],
                                 [-L / 2, H / 2, W],
                                 [L / 2, H / 2, W],
                                 [L / 2, -H / 2, W]])

        for k in range(0, 8):
            plane.AddPart(IGES_Entity116(farfield_pts[k][0], farfield_pts[k][1], farfield_pts[k][2]))

        for k in range(0, 4):
            addLine(plane, farfield_pts, k, (k + 1) % 4)
            addLine(plane, farfield_pts, k + 4, (k + 5) % 4 + 4)
            addLine(plane, farfield_pts, k, k + 4)

        fileout = plane.Generate()
        return fileout


def addLine(container: IGES_Model, pts, i: int, j: int):
    container.AddPart(IGES_Entity110([pts[i][0], pts[i][1], pts[i][2]],
                                     [pts[j][0], pts[j][1], pts[j][2]]))
