import numpy as np
from scipy.interpolate import BSpline
from src.nurbs.utility import equal, to_homogeneous, to_cartesian
from src.nurbs.curve import calc_pnt_param, calc_knot_vector, calc_ctrl_pts
from src.iges.iges_entity128 import IGES_Entity128


class NURBS_Surface(object):
    def __init__(self, U, V, Pw):
        """
        NURBS曲面
        :param U: u方向节点矢量, n+1个元素
        :param V: v方向节点矢量，m+1个元素
        :param Pw: 齐次坐标序列，(n+1)x(m+1)个元素
        """

        self.U = np.copy(U)
        self.V = np.copy(V)
        self.Pw = np.copy(Pw)

        self.n, self.m, dim = Pw.shape
        self.n -= 1
        self.m -= 1

        self.p = len(U) - 1 - self.n - 1
        self.q = len(V) - 1 - self.m - 1

        self.weight = np.zeros((self.n + 1, self.m + 1))
        self.cpt = np.zeros((self.n + 1, self.m + 1, 3))
        self.isPoly = True

        for i in range(0, self.n + 1):
            for j in range(0, self.m + 1):
                self.weight[i][j] = self.Pw[i][j][3]
                if self.isPoly and (not equal(self.weight[i][j], 1.0)):
                    self.isPoly = False
                for k in range(0, 3):
                    self.cpt[i][j][k] = self.Pw[i][j][k] / self.weight[i][j]

    def __call__(self, u, v, k=0, l=0):
        R = []
        for i in range(0, self.n):
            spl = BSpline(self.V, self.Pw[i], self.q)
            R.append(spl(v, l))
        Rw = np.copy(R)
        spl = BSpline(self.U, Rw, self.p)
        pw = spl(u, self.p)
        return to_cartesian(pw)

    def to_iges(self, closed_u, closed_v, periodic_u, periodic_v, form=0):
        return IGES_Entity128(self.U, self.V, self.p, self.q, self.n, self.m, self.cpt, self.weight,
                              closed_u, closed_v, (1 if self.isPoly else 0), periodic_u, periodic_v,
                              self.U[0], self.U[-1], self.V[0], self.V[-1], form)


class GlobalInterpolatedSurf(NURBS_Surface):
    def __init__(self, pts, p, q, umethod='centripetal', vmethod='chord'):
        """
        (n+1)x(m+1)个数据点全局插值，非有理
        不能很好处理局部数据点共面，需小心使用
        :param pts: 待插值数据点
        :param p: u方向次数
        :param q: v方向次数
        :param umethod: u方向参数计算方法
        :param vmethod: v方向参数计算方法 
        """

        n, m, dim = pts.shape
        n -= 1
        m -= 1

        U = np.zeros(n + 1)
        V = np.zeros(m + 1)
        U[-1] = 1.0
        V[-1] = 1.0

        '''Parameters of U direction'''
        dist = np.zeros((n + 1, m + 1))
        for j in range(0, m + 1):
            td = calc_pnt_param(pts[:, j], umethod)
            for i in range(0, n + 1):
                dist[i][j] = td[i]
        for i in range(0, n):
            U[i] = np.mean(dist[i])

        '''Parameters of V Direction'''
        for i in range(0, n + 1):
            td = calc_pnt_param(pts[i], vmethod)
            for j in range(0, m + 1):
                dist[i][j] = td[j]
        for j in range(0, m):
            V[j] = np.mean(dist[:, j])

        '''Knot Vectors'''
        uknot = calc_knot_vector(U, p)
        vknot = calc_knot_vector(V, q)

        '''Control Points'''
        R = np.zeros((n + 1, m + 1, dim))
        for j in range(0, m + 1):
            tp = calc_ctrl_pts(uknot, p, pts[:, j], U)
            for i in range(0, n + 1):
                R[i][j] = tp[i]

        P = np.zeros((n + 1, m + 1, dim))
        for i in range(0, n + 1):
            P[i] = calc_ctrl_pts(vknot, q, R[i], V)

        Pw = np.zeros((n + 1, m + 1, dim + 1))
        for i in range(0, n + 1):
            for j in range(0, m + 1):
                Pw[i][j] = to_homogeneous(P[i][j])

        super(GlobalInterpolatedSurf, self).__init__(uknot, vknot, Pw)


class BilinearSurf(NURBS_Surface):
    def __init__(self, P):
        """
        双线性曲面
        :param P:4个角点, 2x2
        """

        U = np.array([0, 0, 1, 1], float)
        V = np.array([0, 0, 1, 1], float)

        ul, vl, dim = P.shape
        assert ul == 2 and vl == 2

        Pw = np.ones((ul, vl, 4), float)
        for i in range(ul):
            for j in range(vl):
                for d in range(dim):
                    Pw[i][j][d] = P[i][j][d]

        super(BilinearSurf, self).__init__(U, V, Pw)

    def to_iges(self):
        return IGES_Entity128(self.U, self.V, self.p, self.q, self.n, self.m, self.cpt, self.weight,
                              0, 0, (1 if self.isPoly else 0), 0, 0, self.U[0], self.U[-1], self.V[0], self.V[-1], 0)


class RuledSurf(NURBS_Surface):
    def __init__(self):
        pass


class Coons(NURBS_Surface):
    def __init__(self):
        pass


class Skinning(NURBS_Surface):
    def __init__(self):
        pass
