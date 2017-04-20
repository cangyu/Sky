import math
import numpy as np
from scipy.misc import comb
from src.iges.iges_entity128 import IGES_Entity128


def DegreeElevateCurve(n, p, U, Pw, t, MAX_TMP=3000):
    """
    将曲线的次数升高t次
    :param n: Index of the last element in Pw
    :param p: Order of original NURBS Curve
    :param U: Original Knots
    :param Pw: Original Control Pts
    :param t: level to be promoted
    :return: Description of elevated NURBS curve
    """

    m = n + p + 1
    ph = p + t
    ph2 = ph / 2
    dim = Pw[0].shape[0]

    # 新控制点的最后一个下标
    nh = 0

    # 新的节点矢量
    Uh = np.zeros(MAX_TMP, float)

    # 新的控制点
    Qw = np.zeros((MAX_TMP, dim), float)

    bezalfs = np.zeros((p + t + 1, p + 1), float)
    bpts = np.zeros(p + 1, float)
    ebpts = np.zeros(p + t + 1, float)
    Nextbpts = np.zeros(p - 1, float)
    alfs = np.zeros(p - 1, float)

    # 计算升阶的次数
    bezalfs[0][0] = bezalfs[ph][p] = 1.0

    for i in range(1, ph2 + 1):
        inv = 1.0 / comb(ph, i)
        mpi = np.amin([p, i])
        for j in range(np.amax([0, i - t]), mpi + 1):
            bezalfs[i][j] = inv * comb(p, j) * comb(t, i - j)

    for i in range(ph2 + 1, ph):
        mpi = np.amin([p, i])
        for j in range(np.amax([0, i - t]), mpi + 1):
            bezalfs[i][j] = bezalfs[ph - i][p - j]

    mh = ph
    kind = ph + 1
    r = -1
    a = p
    b = p + 1
    cind = 1
    ua = U[0]
    Qw[0] = Pw[0]
    for i in range(0, ph + 1):
        Uh[i] = ua

    # 初始化第一个bezier段
    for i in range(0, p + 1):
        bpts[i] = Pw[i]

    # 沿节点矢量进行循环
    while b < m:
        i = b
        while b < m and U[b] == U[b + 1]:
            b += 1
        mul = b - i + 1
        mh += mul + t
        ub = U[b]
        oldr = r
        r = p - mul

        # 插入节点u(b) r次
        if oldr > 0:
            lbz = (oldr + 2) / 2
        else:
            lbz = 1

        if r > 0:
            rbz = ph - (r + 1) / 2
        else:
            rbz = ph

        if r > 0:
            # 插入节点以获得bezier曲线段
            numer = ub - ua
            for k in range(p, mul, -1):
                alfs[k - mul - 1] = numer / (U[a + k] - ua)
            for j in range(1, r + 1):
                save = r - j
                s = mul + j
                for k in range(p, s - 1, -1):
                    bpts[k] = alfs[k - s] * bpts[k] + (1.0 - alfs[k - s]) * bpts[k - 1]
                Nextbpts[save] = bpts[p]

        # 对bezier曲线段升阶
        for i in range(lbz, ph + 1):
            ebpts[i] = 0.0
            mpi = np.amin([p, i])
            for j in range(np.amax([0, i - t]), mpi + 1):
                ebpts[i] += bezalfs[i][j] * bpts[j]

        if oldr > 1:
            first = kind - 2
            last = kind
            den = ub - ua
            bet = (ub - Uh[kind - 1]) / den
            for tr in range(1, oldr):
                i = first
                j = last
                kj = j - kind + 1
                while j - i > tr:
                    if i < cind:
                        alf = (ub - Uh[i]) / (ua - Uh[i])
                        Qw[i] = alf * Qw[i] + (1.0 - alf) * Qw[i - 1]
                    if j >= lbz:
                        if j - tr <= kind - ph + oldr:
                            gam = (ub - Uh[j - tr]) / den
                            ebpts[kj] = gam * ebpts[kj] + (1.0 - gam) * ebpts[kj + 1]
                        else:
                            ebpts[kj] = bet * ebpts[kj] + (1.0 - bet) * ebpts[kj + 1]
                    i += 1
                    j -= 1
                    kj -= 1
                first -= 1
                last += 1

        # 消去节点u=U[a]结束
        if a != p:
            for i in range(0, ph - oldr):
                Uh[kind] = ua
                kind += 1

        # 将控制点存入Qw
        for j in range(lbz, rbz + 1):
            Qw[cind] = ebpts[j]
            cind += 1

        if b < m:
            # 为下一次循环做准备
            for j in range(0, r):
                bpts[j] = Nextbpts[j]
            for j in range(r, p + 1):
                bpts[j] = Pw[b - p + j]
            a = b
            b += 1
            ua = ub
        else:
            for i in range(0, ph + 1):
                Uh[kind + i] = ub

    nh = mh - ph - 1
    knots = np.zeros(mh + 1, float)
    ctrl_pts = np.zeros((nh + 1, dim), float)
    for i in range(0, mh + 1):
        knots[i] = Uh[i]
    for i in range(0, nh + 1):
        for j in range(0, dim):
            ctrl_pts[i][j] = Qw[i][j]

    return nh, knots, ctrl_pts


def merge_knot(lhs, rhs):
    i = 0
    j = 0
    ret = []
    while i != len(lhs) and j != len(rhs):
        if math.fabs(lhs[i] - rhs[j]) > 1e-7:
            ret.append(lhs[i])
            ret.append(rhs[j])
        else:
            ret.append(lhs[i])

        i += 1
        j += 1

    while i != len(lhs):
        ret.append(lhs[i])
        i += 1

    while j != len(rhs):
        ret.append((rhs[j]))
        j += 1

    return ret


def merge_all_knot(knot_list):
    assert len(knot_list) != 0

    ccu = knot_list[0]
    for i in range(1, len(knot_list)):
        ccu = merge_knot(ccu, knot_list[i])

    return ccu


def GetDistance(lhs, rhs):
    dim = len(lhs)
    assert dim == len(rhs)

    l = 0
    for i in range(0, dim):
        l += math.pow(lhs[i] - rhs[i], 2)

    return math.sqrt(l)


def FindSpan(n, p, u, U):
    """
    确定参数u所在的节点区间的下标
    :param n: 控制点个数=n+1
    :param p: NURBS曲线次数
    :param u: 待考察参数
    :param U: Knots
    :return: u所在的区间下标
    """
    if u == U[n + 1]:  # Corner case: 对于$u=u_m$这一特殊情形，将其节点区间的下标设为n
        return n

    low = p
    high = n + 1
    mid = int((high + low) / 2)
    while u < U[mid] or u >= U[mid + 1]:
        if u < U[mid]:
            high = mid
        else:
            low = mid
        mid = (low + high) / 2

    return mid


def RefineKnotVectCurve(n, p, U, Pw, X, r, Ubar, Qw):
    """
    细化曲线的节点矢量
    :param n: Index of the last element in Pw
    :param p: Order of original NURBS Curve
    :param U: Original Knots:${U_0, U_1, ... , U_m}$
    :param Pw: Original Control Pts:${P_0, P_1, ... , P_n}$
    :param X: Knots to be inserted:${X_0, X_1, ... , X_r}$
    :param r: Index of the last element in X
    :param Ubar: New Knots
    :param Qw: New Control Pts
    :return: None
    """
    if r < 0:
        for i in range(0, len(U)):
            Ubar[i] = U[i]
        for i in range(0, len(Pw)):
            for j in range(0, 3):
                Qw[i][j] = Pw[i][j]
        return

    m = n + p + 1
    a = FindSpan(n, p, X[0], U)
    b = FindSpan(n, p, X[r], U)
    b += 1

    for j in range(0, a - p + 1):
        Qw[j] = Pw[j]

    for j in range(b - 1, n + 1):
        Qw[j + r + 1] = Pw[j]

    for j in range(0, a + 1):
        Ubar[j] = U[j]

    for j in range(b + p, m + 1):
        Ubar[j + r + 1] = U[j]

    i = b + p - 1
    k = b + p + r

    j = r
    while j >= 0:
        while X[j] <= U[i] and i > a:
            Qw[k - p - 1] = Pw[i - p - 1]
            Ubar[k] = U[i]
            k -= 1
            i -= 1

        Qw[k - p - 1] = Qw[k - p]

        for l in range(1, p + 1):
            ind = k - p + 1
            alfa = Ubar[k + 1] - X[j]
            if math.fabs(alfa) == 0.0:
                Qw[ind - 1] = Qw[ind]
            else:
                alfa /= (Ubar[k + 1] - U[i - p + l])
                Qw[ind - 1] = alfa * Qw[ind - 1] + (1.0 - alfa) * Qw[ind]

        Ubar[k] = X[j]
        k -= 1
        j -= 1


def get_diff_knots(lhs, rhs):
    """
    类似集合的差运算
    :param lhs: 装有节点向量的列表，要已经有序！
    :param rhs: 装有节点向量的列表，要已经有序！
    :return: 装有节点向量的列表，其中节点在lhs中，但不在rhs中
    """
    ans = []
    i = 0
    j = 0
    n = len(lhs)
    m = len(rhs)

    while i < n and j < m:
        if math.fabs(lhs[i] - rhs[j]) < 1e-7:
            i += 1
            j += 1
        else:
            ans.append(lhs[i])
            i += 1

    while i < n:
        ans.append(lhs[i])
        i += 1

    return ans


class Skining(object):
    """
    蒙面
    """

    def __init__(self, _crvs):
        self.num = len(_crvs)
        self.crv = _crvs
        self.new_crv_list = []
        self.unified_crv_list = []
        self.u = []
        self.v = []

    def promote(self, p):
        p_list = []
        for crv in self.crv:
            p_list.append(crv.p)

        cpm = np.amax(p_list)
        assert cpm <= p

        self.new_crv_list = []
        for crv in self.crv:
            if p != crv.p:
                self.new_crv_list.append(DegreeElevateCurve(crv.n, crv.p, crv.knots, crv.ctrl_pts, p - crv.p))
            else:
                self.new_crv_list.append([crv.n, crv.knots, crv.ctrl_pts])

    def calc_u_knots(self):
        kl = []
        for i in range(0, self.num):
            kl.append(self.new_crv_list[i][1])

        self.u = merge_all_knot(kl)

    def unify_all_crv(self):
        m = len(self.u) - 1
        self.unified_crv_list = []

        for i in range(0, self.num):
            n, knots, ctrl_pts = self.new_crv_list[i]
            dim = ctrl_pts[0].shape[0]
            diff_knots = get_diff_knots(self.u, knots)
            dfkc = len(diff_knots) - 1
            p = m - n - 1
            new_knots = np.zeros(m + 1, float)
            new_ctrl_pts = np.zeros((n + 1, dim), float)
            RefineKnotVectCurve(len(ctrl_pts) - 1, p, knots, ctrl_pts, diff_knots, dfkc, new_knots, new_ctrl_pts)
            self.unified_crv_list.append((new_knots, new_ctrl_pts))

    def calc_v_knots(self, q):
        K = self.num
        N = len(self.unified_crv_list[0][1])
        d = np.zeros(N, float)
        seg_len = np.zeros((N, K - 1), float)

        for i in range(0, N):
            for j in range(1, K):
                seg_len[i][j - 1] = GetDistance(self.unified_crv_list[j][1][i], self.unified_crv_list[j - 1][1][i])
                d[i] += seg_len[i][j - 1]

        v = np.zeros(K, float)
        v[K - 1] = 1.0
        for k in range(0, K - 2):
            tmp = 0.0
            for i in range(0, N):
                tmp += seg_len[i][k] / d[i]
            v[k + 1] = v[k] + tmp / N

        m = K + q
        knots = np.zeros(m + 1, float)
        for i in range(0, q + 1):
            knots[i] = 0.0
            knots[m - i] = 1.0

        tmp = 0.0
        for i in range(0, q):
            tmp += v[i]

        for j in range(1, K - q):
            tmp -= v[j - 1]
            tmp += v[j + q - 1]
            knots[j + q] = tmp / q

        self.v = knots

    def generate(self, p=3, q=3):
        self.promote(p)
        self.calc_u_knots()
        self.unify_all_crv()
        self.calc_v_knots(q)

        n1 = len(self.u) - 2 - p
        n2 = len(self.v) - 2 - q

        cpts = np.zeros((n1 + 1, n2 + 1, 3), float)
        ws = np.ones((n1 + 1, n2 + 1), float)

        for j in range(0, n2 + 1):
            for i in range(0, n1 + 1):
                for k in range(0, 3):
                    cpts[i][j][k] = self.unified_crv_list[j][1][i][k]

        return IGES_Entity128(self.u, self.v, p, q, n1, n2, cpts, ws)


def LocalSurfInterp(n: int, m: int, Q):
    """
    双三次局部曲面插值
    :param n: u方向最后一个数据点下标
    :param m: v方向最后一个数据点下标
    :param Q: 数据点
    :return: NURBS representation of interpolated surface
    """
    td = np.zeros((n + 1, m + 1, 3), float)
    ub = np.zeros(n + 1, float)
    vb = np.zeros(m + 1, float)
    r = np.zeros(m + 1, float)
    s = np.zeros(n + 1, float)
    U = np.zeros(2 * n + 6, float)
    V = np.zeros(2 * m + 6, float)
    P = np.zeros((2 * n + 2, 2 * m + 2, 3), float)

    total = 0.0
    for l in range(0, m + 1):
        for k in range(1, n + 1):
            d = GetDistance(Q[k][l], Q[k - 1][l])
            ub[k] += d
            r[l] += d
        total += r[l]

    for k in range(1, n):
        ub[k] = ub[k - 1] + ub[k] / total

    ub[n] = 1.0

    total = 0.0
    for k in range(0, n + 1):
        for l in range(1, m + 1):
            d = GetDistance(Q[k][l], Q[k][l - 1])
            vb[l] += d
            s[k] += d
        total += s[k]

    for l in range(1, m):
        vb[l] = vb[l - 1] + vb[l] / total

    vb[m] = 1.0

    # 载入节点矢量U,V
    for i in range(0, 4):
        U[-1 - i] = 1.0
        V[-1 - i] = 1.0

    ci = 4
    for i in range(1, n):
        U[ci + 1] = U[ci] = ub[i]
        ci += 2

    cj = 4
    for j in range(1, m):
        V[cj + 1] = V[ci] = vb[j]
        cj += 2

    return U, V, P  # u方向节点矢量, v方向节点矢量, 控制点


class GordonSurface(object):
    def __init__(self, crv_u_list, crv_v_list):
        pass
