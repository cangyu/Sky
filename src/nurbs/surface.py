import numpy as np
from numpy.linalg import norm
from scipy.interpolate import BSpline
from src.nurbs.utility import equal, to_homogeneous, to_cartesian
from src.nurbs.curve import calc_pnt_param, calc_knot_vector, calc_ctrl_pts, NURBS_Curve
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

    @property
    def n(self):
        """
        U方向最后一个控制点下标
        """

        return self.Pw.shape[0] - 1

    @property
    def m(self):
        """
        V方向最后一个控制点下标
        """

        return self.Pw.shape[1] - 1

    @property
    def p(self):
        """
        U方向次数
        """

        return len(self.U) - self.n - 2

    @property
    def q(self):
        """
        V方向次数
        """

        return len(self.V) - self.m - 2

    def __call__(self, u, v, k=0, l=0, return_cartesian=True):
        R = []
        for i in range(0, self.n + 1):
            spl = BSpline(self.V, self.Pw[i], self.q)
            R.append(spl(v, l))
        Rw = np.copy(R)
        spl = BSpline(self.U, Rw, self.p)
        pw = spl(u, k)
        return to_cartesian(pw) if return_cartesian else pw

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


class ExtrudedSurf(NURBS_Surface):
    def __init__(self, crv: NURBS_Curve, dir):
        """
        拉伸曲面
        :param crv: Curve to be extruded.
        :param dir: Direction vector.
        """

        U = np.copy(crv.U)
        V = np.array([0, 0, 1, 1], float)
        n = len(crv.cpt)
        Pw = np.zeros((n, 2, 4))
        for i in range(n):
            Pw[i][0] = Pw[i][1] = np.copy(crv.Pw[i])
            wdir = to_homogeneous(dir, Pw[i][0][3])
            for d in range(3):
                Pw[i][1][d] += wdir[d]

        super(ExtrudedSurf, self).__init__(U, V, Pw)


class RuledSurf(NURBS_Surface):
    def __init__(self, c1, c2):
        """
        生成V方向的直纹面,即两条曲线之间的线性插值
        :param c1: 第1条曲线
        :type c1: NURBS_Curve
        :param c2: 第2条曲线
        :type c2:NURBS_Curve
        """

        '''Check'''
        if not equal(c1.U[0], c2.U[0]):
            raise ValueError('Incompatible starting knot!')
        if not equal(c1.U[-1], c2.U[-1]):
            raise ValueError('Incompatible ending knot!')

        '''Knot vector'''
        p = max(c1.p, c2.p)
        c1.elevate(p - c1.p, self_update=True, return_raw=True)
        c2.elevate(p - c2.p, self_update=True, return_raw=True)

        if not equal(norm(c1.U - c2.U), 0):
            all_knot = merge_knot(c1.U, c2.U)
            x1 = different_knot(c1.U, all_knot)
            x2 = different_knot(c2.U, all_knot)
            c1.refine(x1)
            c2.refine(x2)

        uknot = c1.U
        vknot = np.array([0, 0, 1, 1], float)

        '''Control points'''
        pw = np.zeros((len(c1.Pw), 2, 4))
        for i in range(len(c1.Pw)):
            pw[i][0] = np.copy(c1.Pw[i])
            pw[i][1] = np.copy(c2.Pw[i])

        super(RuledSurf, self).__init__(uknot, vknot, pw)


class Coons(NURBS_Surface):
    def __init__(self):
        pass


class Skinning(NURBS_Surface):
    def __init__(self, crv, p, q, vmethod='chord'):
        """
        蒙皮曲面，非有理
        :param crv: 非有理曲线集合
        :param p: 目标曲面u方向次数(曲线方向)
        :param q: 目标曲面v方向次数(展向)
        :param vmethod: v方向插值方法
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
    ans = np.zeros(int(sum(cnt)))
    k = 0
    for i in range(0, len(val)):
        for j in range(0, cnt[i]):
            ans[k] = val[i]
            k += 1

    return ans
