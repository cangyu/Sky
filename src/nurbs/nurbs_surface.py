import numpy as np
import math
from scipy.misc import comb
from src.nurbs.basis import *
from src.nurbs.nurbs_curve import *
from src.iges.iges_entity128 import IGES_Entity128


class NURBS_Surface(object):
    def __init__(self, U, V, Pw):
        """
        NURBS曲面
        :param U: u方向节点矢量
        :param V: v方向节点矢量
        :param Pw: 齐次坐标序列
        """

        self.U = np.copy(U)
        self.V = np.copy(V)
        self.Pw = np.copy(Pw)

        self.n, self.m, dim = Pw.shape
        self.n -= 1
        self.m -= 1

        self.p = len(U) - 1 - self.n - 1
        self.q = len(V) - 1 - self.m - 1

        self.N1 = Basis(U, self.p)
        self.N2 = Basis(V, self.q)

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

    def w(self, u, v, k=0, l=0):
        ans = 0.0
        nipd = np.zeros(self.n + 1)
        mjqd = np.zeros(self.m + 1)
        for i in range(0, self.n + 1):
            nipd[i] = self.N1(i, self.p, u, k)
        for j in range(0, self.m + 1):
            mjqd[j] = self.N2(j, self.q, v, l)

        for i in range(0, self.n + 1):
            for j in range(0, self.m + 1):
                ans += nipd[i] * mjqd[j] * self.Pw[i][j][3]

        return ans

    def __call__(self, u, v, k=0, l=0):
        nipd = np.zeros(self.n + 1)
        mjqd = np.zeros(self.m + 1)
        for i in range(0, self.n + 1):
            nipd[i] = self.N1(i, self.p, u, k)
        for j in range(0, self.m + 1):
            mjqd[j] = self.N2(j, self.q, v, l)

        Akl = np.zeros(3)
        for i in range(0, self.n + 1):
            for j in range(0, self.m + 1):
                for kk in range(0, 3):
                    Akl[kk] += nipd[i] * mjqd[j] * self.Pw[i][j][kk]

        ccw1 = np.zeros(3)
        i = 1
        while i <= k:
            ccw1 += comb(k, i, exact=True) * self.w(u, v, i, 0) * self.__call__(u, v, k - i, l)
            i += 1

        ccw2 = np.zeros(3)
        j = 1
        while j <= l:
            ccw2 += comb(l, j, exact=True) * self.w(u, v, 0, j) * self.__call__(u, v, k, l - j)
            j += 1

        ccw3 = np.zeros(3)
        i = 1
        while i <= k:
            ci = comb(k, i, exact=True)
            j = 1
            while j <= l:
                ccw3 += ci * comb(l, j, exact=True) * self.w(u, v, i, j) * self.__call__(u, v, k - i, l - j)
                j += 1
            i += 1

        wuv = self.w(u, v)
        if equal(wuv, 0.0):
            raise ValueError("Invalid weight settings!")

        return (Akl - ccw1 - ccw2 - ccw3) / wuv

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


class Skinning(NURBS_Surface):
    def __init__(self, crv, p, q, vmethod='chord'):
        """
        蒙皮曲面，非有理
        :param crv: 非有理曲线集合
        :param p: 目标曲面u方向次数(曲线方向)
        :param q: 目标曲面v方向次数(展向)
        :param vmethod: v方向插值方法
        :return: None
        """

        '''Promote all curve to p order'''
        self.crv = []
        for i in range(0, len(crv)):
            self.crv.append(crv[i])

        for i in range(0, len(self.crv)):
            self.crv[i].elevate(p - self.crv[i].p)

        '''Merge all knot vectors in U direction'''
        uknot = np.copy(self.crv[0].U)
        for i in range(1, len(self.crv)):
            ck = np.copy(self.crv[i].U)
            uknot = merge_knot(uknot, ck)

        '''Unify all curve knot vector'''
        for i in range(0, len(self.crv)):
            X = different_knot(uknot, self.crv[i].U)
            for u in X:
                self.crv[i].insert_knot(u)

        '''Knot vector in V direction'''
        n = len(uknot) - 1 - p - 1
        m = len(self.crv) - 1
        pnt = np.zeros((n + 1, m + 1, 3))
        for j in range(0, m + 1):
            for i in range(0, n + 1):
                pnt[i][j] = to_cartesian(self.crv[j].Pw[i])

        vparam = np.zeros((n + 1, m + 1))
        vk = []
        for i in range(0, n + 1):
            vparam[i] = calc_pnt_param(pnt[i], vmethod)
            vk.append(calc_knot_vector(vparam[i], q))

        vknot = np.mean(vk, axis=0)

        '''Calculate control points'''
        Q = np.zeros((n + 1, m + 1, 3))
        Qw = np.zeros((n + 1, m + 1, 4))

        for i in range(0, n + 1):
            Q[i] = calc_ctrl_pts(vknot, q, pnt[i], vparam[i])

        for i in range(0, n + 1):
            for j in range(0, m + 1):
                Qw[i][j] = to_homogeneous(Q[i][j])

        super(Skinning, self).__init__(uknot, vknot, Qw)


def merge_knot(lhs, rhs):
    """
    合并两个节点矢量
    :param lhs: 第1个节点矢量
    :param rhs: 第2个节点矢量
    :return: 合并后的节点矢量, lhs union rhs
    """

    lval, lcnt = np.unique(lhs, return_counts=True)
    rval, rcnt = np.unique(rhs, return_counts=True)

    val = np.concatenate((lval, rval))
    val = np.unique(val)
    ans = []

    for v in val:
        if v in lval and v in rval:
            li = np.searchsorted(lval, v)
            ri = np.searchsorted(rval, v)
            cc = max(lcnt[li], rcnt[ri])
            for i in range(0, cc):
                ans.append(v)
        else:
            if v in lval:
                li = np.searchsorted(lval, v)
                for i in range(0, lcnt[li]):
                    ans.append(v)
            else:
                ri = np.searchsorted(rval, v)
                for i in range(0, rcnt[ri]):
                    ans.append(v)

    ret = np.zeros(len(ans), float)
    for i in range(0, len(ret)):
        ret[i] = ans[i]

    return ret


def different_knot(lhs, rhs):
    """
    求两个节点矢量中的不同部分
    :param lhs: 第1个节点矢量
    :param rhs: 第2个节点矢量
    :return: lhs subtract rhs
    """

    lval, lcnt = np.unique(lhs, return_counts=True)
    rval, rcnt = np.unique(rhs, return_counts=True)

    '''Count each'''
    val = []
    cnt = []
    for i in range(0, len(lval)):
        if lval[i] in rval:
            k = np.searchsorted(rval, lval[i])
            lvc = lcnt[i]
            rvc = rcnt[k]
            if lvc > rvc:
                val.append(lval[i])
                cnt.append(lvc - rvc)
        else:
            val.append(lval[i])
            cnt.append(lcnt[i])

    '''Assemble'''
    ans = np.zeros(int(np.sum(cnt)))
    k = 0
    for i in range(0, len(val)):
        for j in range(0, cnt[i]):
            ans[k] = val[i]
            k += 1

    return ans
