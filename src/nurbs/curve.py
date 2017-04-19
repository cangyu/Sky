import numpy as np
import scipy.linalg
from scipy import interpolate
from scipy.special import factorial
import math
from src.iges.iges_entity112 import IGES_Entity112


def BasisFuns(i, u, p, U):
    """
    计算在给定点u，所有p次非零B样条基函数的值(p+1个):$N_{i,p}(u)$
    （只有 $N_{i-p,p}(u) -- N_{i,p}(u)$ 不为0）
    """
    left = np.zeros(p + 1, float)
    right = np.zeros(p + 1, float)
    N = np.zeros(p + 1, float)

    N[0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - U[i + 1 - j]
        right[j] = U[i + j] - u
        saved = 0.0
        for r in range(0, j):
            tmp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * tmp
            saved = left[j - r] * tmp

        N[j] = saved

    return N


class Curve(object):
    def __init__(self, pts, p=5):
        """
        p阶，n+1个控制点，m+1个节点,全局插值，非有理
        """
        cur_shape = pts[0].shape
        for i in range(1, len(pts)):
            if cur_shape != pts[i].shape:
                raise Exception("Inconsistent shape!")

        self.dim = cur_shape[0]
        self.p = p
        self.n = len(pts) - 1
        self.m = self.n + self.p + 1

        self.pts = np.zeros((self.n + 1, self.dim), float)
        for i in range(0, self.n + 1):
            for j in range(0, self.dim):
                self.pts[i][j] = pts[i][j]

        self.param = np.zeros(self.n + 1, float)
        self.knots = np.zeros(self.m + 1, float)
        self.coef = np.zeros((self.n + 1, self.n + 1), float)
        self.ctrl_pts = np.zeros((self.n + 1, self.dim), float)
        self.weights = np.ones(self.n + 1, float)

    def calc_param(self, method='centripetal'):
        self.param[0] = 0.0
        self.param[self.n] = 1.0

        if method not in ['uniform', 'chord', 'centripetal']:
            raise Exception("Invalid parameter!")

        if method == 'uniform':
            inc = 1 / self.n
            for i in range(1, self.n):
                self.param[i] = self.param[i - 1] + inc
        else:
            dist = np.zeros(self.n + 1, float)
            for i in range(1, self.n + 1):
                for j in range(0, self.dim):
                    dist[i] += math.pow(self.pts[i][j] - self.pts[i - 1][j], 2)
                dist[i] = math.sqrt(dist[i])

            d = 0
            if method == 'chord':
                for i in range(1, self.n + 1):
                    d += dist[i]
            else:
                for i in range(1, self.n + 1):
                    dist[i] = math.sqrt(dist[i])
                    d += dist[i]

            for i in range(1, self.n):
                self.param[i] = self.param[i - 1] + dist[i] / d

    def calc_knots(self):
        """
        取平均值方法计算节点
        """
        for i in range(0, self.p + 1):
            self.knots[i] = 0.0
            self.knots[self.m - i] = 1.0

        acc = 0.0
        for i in range(0, self.p):
            acc += self.param[i]

        for j in range(1, self.n - self.p + 1):
            acc -= self.param[j - 1]
            acc += self.param[self.p - 1 + j]
            self.knots[self.p + j] = acc / self.p

    def calc_coef(self):
        """
        计算系数矩阵
        """
        rg = np.zeros(self.n + 1, int)
        rg[0] = self.p
        rg[self.n] = self.n

        for i in range(1, self.n):
            left = self.p
            right = self.n

            while left < right:
                mid = int((left + right + 1) / 2)
                if self.param[i] < self.knots[mid]:
                    right = mid - 1
                else:
                    left = mid

            rg[i] = left

        for k in range(0, self.n + 1):
            i = rg[k]
            u = self.param[k]

            N = BasisFuns(i, u, self.p, self.knots)

            for j in range(i - self.p, i + 1):
                self.coef[k][j] = N[j - i + self.p]

    def calc_ctrl_pts(self):
        """
        解self.dim个线性方程组反求得控制点坐标，每个方程组中有self.n+1个方程
        """
        Q = np.zeros((self.dim, self.n + 1), float)
        P = np.zeros((self.dim, self.n + 1), float)

        for i in range(0, self.dim):
            for j in range(0, self.n + 1):
                Q[i][j] = self.pts[j][i]

        for i in range(0, self.dim):
            P[i] = scipy.linalg.solve(self.coef, Q[i])

        for i in range(0, self.n + 1):
            for j in range(0, self.dim):
                self.ctrl_pts[i][j] = P[j][i]

    def generate(self, method='centripetal'):
        self.calc_param(method)
        self.calc_knots()
        self.calc_coef()
        self.calc_ctrl_pts()

        # 节点矢量， 权系数， 控制点
        return self.knots, self.weights, self.ctrl_pts


class Spline(object):
    def __init__(self, _u, _x, _y, _z, _order=3):

        self.order = _order
        self.u = np.zeros(len(_u), float)
        self.x = np.zeros(len(_x), float)
        self.y = np.zeros(len(_y), float)
        self.z = np.zeros(len(_z), float)

        for i in range(0, len(_u)):
            self.u[i] = _u[i]

        for i in range(0, len(_x)):
            self.x[i] = _x[i]

        for i in range(0, len(_y)):
            self.y[i] = _y[i]

        for i in range(0, len(_z)):
            self.z[i] = _z[i]

    def generate(self, bc_x=([(2, 0)], [(2, 0)]), bc_y=([(2, 0)], [(2, 0)]), bc_z=([(2, 0)], [(2, 0)])):
        """
        We take Natural BC for each direction by default.
        :param bc_x: BC for X-direction
        :param bc_y: BC for Y-direction
        :param bc_z: BC for Z-direction
        :return: An IGES entity representing parametric curve
        """
        # Coefficient Matrix
        self.cm = np.zeros((len(self.u), 3, 4))

        # Interpolation Function
        fx = interpolate.make_interp_spline(self.u, self.x, k=self.order, bc_type=bc_x)
        fy = interpolate.make_interp_spline(self.u, self.y, k=self.order, bc_type=bc_y)
        fz = interpolate.make_interp_spline(self.u, self.z, k=self.order, bc_type=bc_z)

        f = [fx, fy, fz]
        for k in range(0, len(self.u)):
            for i in range(0, 3):
                for j in range(0, 4):
                    self.cm[k][i][j] = float(f[i](self.u[k], j) / factorial(j))

        return IGES_Entity112(self.u, self.cm)
