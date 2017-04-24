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

    def generate(self, show_airfoil=True, show_surf=False, show_airfoil_cross_lines=True, show_airfoil_pts=False):
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

        fileout = plane.Generate()
        return fileout
