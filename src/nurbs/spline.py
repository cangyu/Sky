import numpy as np
from scipy import interpolate
from scipy.special import factorial
from src.nurbs.utility import pnt_dist
from src.iges.iges_entity112 import IGES_Entity112


class Spline(object):
    def __init__(self):
        """
        3次样条曲线
        """
        self.u = None
        self.cm = None
        self.x = None
        self.y = None
        self.z = None

    def interpolate(self, pts, bc_x=([(2, 0)], [(2, 0)]), bc_y=([(2, 0)], [(2, 0)]), bc_z=([(2, 0)], [(2, 0)])):
        """
        利用scipy中的bspline插值，从给定点序列构造3次样条曲线
        """

        '''Copy parameters and coordinates'''
        n = len(pts)
        self.u = np.zeros(n, float)
        xc = np.zeros(n, float)
        yc = np.zeros(n, float)
        zc = np.zeros(n, float)
        for i in range(0, n):
            xc[i] = pts[i][0]
            yc[i] = pts[i][1]
            zc[i] = pts[i][2]

        '''Natural coordinates'''
        for i in range(1, n):
            self.u[i] = pnt_dist(pts[i], pts[i - 1]) + self.u[i - 1]

        '''Interpolation Function'''
        order = 3
        fx = interpolate.make_interp_spline(self.u, xc, k=order, bc_type=bc_x)
        fy = interpolate.make_interp_spline(self.u, yc, k=order, bc_type=bc_y)
        fz = interpolate.make_interp_spline(self.u, zc, k=order, bc_type=bc_z)

        '''Normalized representation of each dimension'''
        self.x = lambda u: fx(u * self.u[-1])
        self.y = lambda u: fy(u * self.u[-1])
        self.z = lambda u: fz(u * self.u[-1])

        '''Coefficient matrix'''
        f = [fx, fy, fz]
        self.cm = np.zeros((n, 3, order + 1))
        for k in range(0, n):
            for i in range(0, 3):
                for j in range(0, order + 1):
                    self.cm[k][i][j] = f[i](self.u[k], j) / factorial(j)

    def to_iges(self):
        return IGES_Entity112(self.u, self.cm)
