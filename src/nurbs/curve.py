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


def GetDistance(lhs, rhs):
    tmp = 0.0
    for j in range(0, 3):
        tmp += math.pow(lhs[j] - rhs[j], 2)
    return math.sqrt(tmp)


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
    def __init__(self, pts, p=3, smooth=1):
        # Number of pts and Order
        self.n = pts.shape[0]
        self.order = p

        # Copy parameters and coordinates
        self.u = np.zeros(self.n, float)
        self.x = np.zeros(self.n, float)
        self.y = np.zeros(self.n, float)
        self.z = np.zeros(self.n, float)
        for i in range(0, self.n):
            self.x[i] = pts[i][0]
            self.y[i] = pts[i][1]
            self.z[i] = pts[i][2]

        # 用自然坐标参数化
        for i in range(1, self.n):
            self.u[i] = GetDistance(pts[i], pts[i - 1]) + self.u[i - 1]

        # 边界条件
        bc_x = ([(2, 0)], [(2, 0)])
        bc_y = ([(2, 0)], [(2, 0)])
        bc_z = ([(2, 0)], [(2, 0)])

        if smooth:
            pder_x = (self.x[1] - self.x[0]) / self.u[1]
            pder_y = (self.y[1] - self.y[0]) / self.u[1]
            pder_z = (self.z[1] - self.z[0]) / self.u[1]
            bc_x = ([(1, pder_x)], [(1, pder_x)])
            bc_y = ([(1, pder_y)], [(1, pder_y)])
            bc_z = ([(1, pder_z)], [(1, pder_z)])

        # Interpolation Function
        fx = interpolate.make_interp_spline(self.u, self.x, k=self.order, bc_type=bc_x)
        fy = interpolate.make_interp_spline(self.u, self.y, k=self.order, bc_type=bc_y)
        fz = interpolate.make_interp_spline(self.u, self.z, k=self.order, bc_type=bc_z)
        f = [fx, fy, fz]

        # Calculate coefficient matrix
        self.cm = np.zeros((self.n, 3, self.order + 1))
        for k in range(0, self.n):
            for i in range(0, 3):
                for j in range(0, self.order + 1):
                    self.cm[k][i][j] = float(f[i](self.u[k], j) / factorial(j))

    def iges(self):
        assert self.order <= 3
        return IGES_Entity112(self.u, self.cm)
