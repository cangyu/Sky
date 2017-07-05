from numpy.linalg import norm
from scipy.interpolate import BSpline
from src.nurbs.utility import *
from src.nurbs.transform import Quaternion
from src.nurbs.curve import calc_pnt_param, calc_knot_vector, calc_ctrl_pts, ClampedNURBSCrv
from src.iges.iges_entity128 import IGES_Entity128


class ClampedNURBSSurf(object):
    def __init__(self, u, v, pw):
        """
        NURBS曲面
        :param u: u方向节点矢量, n+1个元素
        :param v: v方向节点矢量，m+1个元素
        :param pw: 齐次坐标序列，(n+1)x(m+1)个元素
        """

        self.U = np.copy(u)
        self.V = np.copy(v)
        self.Pw = np.copy(pw)

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

    @property
    def weight(self):
        """
        权系数
        """

        return self.Pw[:, :, -1]

    @property
    def cpt(self):
        """
        不带权控制点
        """

        ans = np.zeros((self.n + 1, self.m + 1, 3))
        for i in range(self.n + 1):
            for j in range(self.m + 1):
                ans[i][j] = to_cartesian(self.Pw[i][j])
        return ans

    def __call__(self, u, v, k=0, l=0, return_cartesian=True):
        """
        求在给定位置(u,v)处的导矢量
        :param u: U方向参数
        :param v: V方向参数
        :param k: U方向求导次数
        :param l: V方向求导次数
        :param return_cartesian: 返回结果形式
        :return: (u,v)处偏导矢量
        """

        r = []
        for i in range(0, self.n + 1):
            spl = BSpline(self.V, self.Pw[i], self.q)
            r.append(spl(v, l))
        rw = np.copy(r)
        spl = BSpline(self.U, rw, self.p)
        pw = spl(u, k)
        return to_cartesian(pw) if return_cartesian else pw

    def reset(self, u, v, pw):
        """
        重置曲面
        :param u: u方向节点矢量, n+1个元素
        :param v: v方向节点矢量，m+1个元素
        :param pw: 齐次坐标序列，(n+1)x(m+1)个元素
        """

        self.U = np.copy(u)
        self.V = np.copy(v)
        self.Pw = np.copy(pw)

    def reverse(self, direction):
        """
        曲面反向
        """

        if direction not in ('U', 'V', 'UV'):
            raise ValueError('Invalid direction choice!')

        if direction in ('U', 'UV'):
            self.U = np.ones(self.U.shape) - self.U
            self.Pw = self.Pw[::-1, :, :]

        if direction in ('V', 'UV'):
            self.V = np.ones(self.V.shape) - self.V
            self.Pw = self.Pw[:, ::-1, :]

    def __repr__(self):
        return 'U Knot:\n{}\nV Knot:\n{}\nControl points:\n{}\n'.format(self.U, self.V, self.Pw)

    def pan(self, delta):
        """
        曲面整体平移
        :param delta: 偏移矢量
        :return: None.
        """

        dv = np.zeros(3)
        array_smart_copy(delta, dv)

        for i in range(self.n + 1):
            for j in range(self.m + 1):
                cv = to_cartesian(self.Pw[i][j]) + dv
                self.Pw[i][j] = to_homogeneous(cv, self.Pw[i][j][-1])

    def rotate(self, ref, ax, ang):
        """
        将曲面绕过指定点的转轴旋转一定角度
        :param ref: 参考点
        :param ax: 旋转轴方向向量，按右手定则确定角度的正方向
        :param ang: 旋转角(Degree)
        :return: None.
        """

        q = Quaternion.from_u_theta(ax, math.radians(ang))
        for i in range(self.n + 1):
            for j in range(self.m + 1):
                cv = to_cartesian(self.Pw[i][j]) - ref
                cv = ref + q.rotate(cv)
                self.Pw[i][j] = to_homogeneous(cv, self.Pw[i][j][-1])

    def to_iges(self, closed_u=0, closed_v=0, periodic_u=0, periodic_v=0, form=0):
        """
        将曲面以IGES标准中第128号实体呈现
        :param closed_u: U方向是否封闭
        :type closed_u: int
        :param closed_v: V方向是否封闭
        :type closed_v: int
        :param periodic_u: U方向是否是周期性的
        :param periodic_v: V方向是否是周期性的
        :param form: IGES中特定形式
        :return: IGES_Entity128 Object
        """

        w = self.weight
        poly = 0 if (w != np.ones(w.shape)).any() else 1
        cpt = self.cpt

        return IGES_Entity128(self.U, self.V, self.p, self.q, self.n, self.m, cpt, w,
                              closed_u, closed_v, poly, periodic_u, periodic_v,
                              self.U[0], self.U[-1], self.V[0], self.V[-1], form)

    def insert_knot(self, uv, r=1, direction='U'):
        """
        曲面插入节点
        :param uv: 待插入节点值
        :param r: 插入次数
        :param direction: 插入的方向
        :return: None
        """

        if direction not in ('U', 'V'):
            raise AssertionError('Invalid direction choice!')

        if direction == 'U':
            crv_list = []
            npw = np.zeros((self.n + 2, self.m + 1, 4))
            for j in range(self.m + 1):
                cc = ClampedNURBSCrv(self.U, self.Pw[:, j])
                cc.insert_knot(uv, r)
                crv_list.append(cc)
                for i in range(self.n + 2):
                    npw[i][j] = np.copy(cc.Pw[i])

            self.reset(crv_list[0].U, self.V, npw)

        else:
            crv_list = []
            npw = np.zeros((self.n + 1, self.m + 2, 4))
            for i in range(self.n + 1):
                cc = ClampedNURBSCrv(self.V, self.Pw[i, :])
                cc.insert_knot(uv, r)
                crv_list.append(cc)
                for j in range(self.m + 2):
                    npw[i][j] = np.copy(cc.Pw[j])

            self.reset(self.U, crv_list[0].V, npw)

    def extract(self, direction, uv):
        """
        提取等参数线
        :param direction: 方向
        :param uv: 等参数值
        :return: 给定方向上的等参数线
        :rtype: ClampedNURBSCrv
        """

        if direction not in ('U', 'V'):
            raise AssertionError('Invalid direction choice!')
        if np.less(uv, 0) or np.greater(uv, 1):
            raise AssertionError('Invalid parameter!')

        if direction == 'U':
            nqw = np.zeros((self.m + 1, 4))
            for j in range(self.m + 1):
                spl = BSpline(self.U, self.Pw[:, j, :], self.p)
                nqw[j] = spl(uv)

            return ClampedNURBSCrv(self.V, nqw)

        else:
            npw = np.zeros(self.n + 1, 4)
            for i in range(self.n + 1):
                spl = BSpline(self.V, self.Pw[i, :, :], self.q)
                npw[i] = spl(uv)

            return ClampedNURBSCrv(self.U, npw)

    def refine(self, direction, extra_knot):
        """
        细化节点矢量
        :param direction: 方向选择
        :param extra_knot: 待插入节点数组
        :return: None
        """

        if direction not in ('U', 'V'):
            raise AssertionError('Invalid direction choice!')
        if len(extra_knot) == 0:
            return

        crv_list = []
        if direction == 'U':
            nh = self.n + 1 + len(extra_knot)
            npw = np.zeros((nh, self.m + 1, 4))
            for j in range(self.m + 1):
                cc = ClampedNURBSCrv(self.U, self.Pw[:, j, :])
                cc.refine(extra_knot)
                crv_list.append(cc)
                for i in range(nh):
                    npw[i][j] = np.copy(cc.Pw[i])

            self.reset(crv_list[0].U, self.V, npw)

        else:
            mh = self.m + 1 + len(extra_knot)
            npw = np.zeros((self.n + 1, mh, 4))
            for i in range(self.n + 1):
                cc = ClampedNURBSCrv(self.V, self.Pw[i, :, :])
                cc.refine(extra_knot)
                crv_list.append(cc)
                for j in range(mh):
                    npw[i][j] = np.copy(cc.Pw[j])

            self.reset(self.U, crv_list[0].U, npw)

    def elevate(self, tu=0, tv=0):
        """
        曲面升阶
        :param tu: U方向升阶次数
        :type tu: int
        :param tv: V方向升阶次数
        :type tv: int
        :return: None.
        """

        if tu < 0 or tv < 0:
            raise AssertionError('Invalid promotion!')

        if tu > 0:
            crv_list = []
            for j in range(self.m + 1):
                cc = ClampedNURBSCrv(self.U, self.Pw[:, j, :])
                cc.elevate(tu)
                crv_list.append(cc)

            nh = len(crv_list[0].Pw)
            npw = np.zeros((nh, self.m + 1, 4))
            for j in range(self.m + 1):
                for i in range(nh):
                    npw[i][j] = np.copy(crv_list[j].Pw[i])

            self.reset(crv_list[0].U, self.V, npw)

        if tv > 0:
            crv_list = []
            for i in range(self.n + 1):
                cc = ClampedNURBSCrv(self.V, self.Pw[i, :, :])
                cc.elevate(tv)
                crv_list.append(cc)

            mh = len(crv_list[0].Pw)
            npw = np.zeros((self.n + 1, mh, 4))
            for i in range(self.n + 1):
                for j in range(mh):
                    npw[i][j] = np.copy(crv_list[i].Pw[j])

            self.reset(self.U, crv_list[0].U, npw)


class GlobalInterpolatedSurf(ClampedNURBSSurf):
    def __init__(self, pts, p, q, u_method='centripetal', v_method='chord'):
        """
        (n+1)x(m+1)个数据点全局插值，非有理
        不能很好处理局部数据点共面，需小心使用
        :param pts: 待插值数据点
        :param p: u方向次数
        :param q: v方向次数
        :param u_method: u方向参数计算方法
        :param v_method: v方向参数计算方法
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
            td = calc_pnt_param(pts[:, j], u_method)
            for i in range(0, n + 1):
                dist[i][j] = td[i]
        for i in range(0, n):
            U[i] = np.mean(dist[i])

        '''Parameters of V Direction'''
        for i in range(0, n + 1):
            td = calc_pnt_param(pts[i], v_method)
            for j in range(0, m + 1):
                dist[i][j] = td[j]
        for j in range(0, m):
            V[j] = np.mean(dist[:, j])

        '''Knot Vectors'''
        u_knot = calc_knot_vector(U, p)
        v_knot = calc_knot_vector(V, q)

        '''Control Points'''
        R = np.zeros((n + 1, m + 1, dim))
        for j in range(0, m + 1):
            tp = calc_ctrl_pts(u_knot, p, pts[:, j], U)
            for i in range(0, n + 1):
                R[i][j] = tp[i]

        P = np.zeros((n + 1, m + 1, dim))
        for i in range(0, n + 1):
            P[i] = calc_ctrl_pts(v_knot, q, R[i], V)

        Pw = np.zeros((n + 1, m + 1, dim + 1))
        for i in range(0, n + 1):
            for j in range(0, m + 1):
                Pw[i][j] = to_homogeneous(P[i][j])

        super(GlobalInterpolatedSurf, self).__init__(u_knot, v_knot, Pw)


class BilinearSurf(ClampedNURBSSurf):
    def __init__(self, P):
        """
        双线性曲面

        ^ V direction
        |
        |P[0][1]        P[1][1]
        ----------------
        |              |
        |              |
        |     SURF     |
        |              |
        |P[0][0]       |P[1][0]
        --------------------------> U direction

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


class ExtrudedSurf(ClampedNURBSSurf):
    def __init__(self, crv, direction):
        """
        拉伸曲面
        :param crv: Curve to be extruded.
        :type crv: ClampedNURBSCrv
        :param direction: Direction vector.
        """

        U = np.copy(crv.U)
        V = np.array([0, 0, 1, 1], float)
        n = len(crv.cpt)
        Pw = np.zeros((n, 2, 4))
        for i in range(n):
            Pw[i][0] = Pw[i][1] = np.copy(crv.Pw[i])
            wdir = to_homogeneous(direction, Pw[i][0][3])
            for d in range(3):
                Pw[i][1][d] += wdir[d]

        super(ExtrudedSurf, self).__init__(U, V, Pw)


class RuledSurf(ClampedNURBSSurf):
    def __init__(self, c1, c2):
        """
        生成V方向的直纹面,即两条曲线之间的线性插值
        :param c1: 第1条曲线
        :type c1: ClampedNURBSCrv
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

        if len(c1.U) != len(c2.U) or not equal(norm(c1.U - c2.U), 0):
            all_knot = merge_knot(c1.U, c2.U)
            x1 = different_knot(all_knot, c1.U)
            x2 = different_knot(all_knot, c2.U)
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


class Coons(ClampedNURBSSurf):
    def __init__(self, c0u, c1u, c0v, c1v):
        """
        双线性混合Coons曲面

         ^ V direction
         |
         |     c1u
         ------->--------
         |              |
         |              |
     c0v ^     SURF     ^ c1v
         |              |
         |              |
         ------->-----------> U direction
               c0u

        :param c0u:沿U方向第1条曲线
        :type c0u: NURBS_Curve
        :param c1u:沿U方向第2条曲线
        :type c1u: NURBS_Curve
        :param c0v:沿V方向第1条曲线
        :type c0v: ClampedNURBSCrv
        :param c1v:沿V方向第2条曲线
        :type c1v: NURBS_Curve
        """

        '''Check 4 corners'''
        assert equal(norm(c0u(0) - c0v(0)), 0.0)
        assert equal(norm(c0u(1) - c1v(0)), 0.0)
        assert equal(norm(c1v(1) - c1u(1)), 0.0)
        assert equal(norm(c0v(1) - c1u(0)), 0.0)

        s = np.zeros((2, 2, 3))
        s[0][0] = np.copy(c0u(0))
        s[0][1] = np.copy(c0v(1))
        s[1][0] = np.copy(c1v(0))
        s[1][1] = np.copy(c1u(1))

        r1 = RuledSurf(c0u, c1u)
        r2 = RuledSurf(c0v, c1v)
        t = BilinearSurf(s)


class Skinning(ClampedNURBSSurf):
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
    val = np.unique(np.concatenate((lval, rval)))

    ans = []
    for v in val:
        if v in lval and v in rval:
            li = np.searchsorted(lval, v)
            ri = np.searchsorted(rval, v)
            cc = max(lcnt[li], rcnt[ri])
            for i in range(cc):
                ans.append(v)
        else:
            if v in lval:
                li = np.searchsorted(lval, v)
                for i in range(lcnt[li]):
                    ans.append(v)
            else:
                ri = np.searchsorted(rval, v)
                for i in range(rcnt[ri]):
                    ans.append(v)

    return np.copy(ans)


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
