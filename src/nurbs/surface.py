import math
import numpy as np
from scipy.misc import comb


def FindSpan(n, p, u, U):
    """
    确定参数u所在的节点区间的下标
    :param n: 控制点个数=n+1
    :param p: NURBS曲线次数
    :param u: 待考察参数
    :param U: Knots
    :return: u所在的区间下标
    """
    if u == U[n + 1]:  # Corner case: 对于u=u_m这一特殊情形，将其节点区间的下标设为n
        return n

    low = p
    high = n + 1
    mid = low + (high - low) / 2
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
    :param n: 控制点个数=n+1
    :param p: NURBS曲线次数
    :param U: Original Knots
    :param Pw: Original Control Pts
    :param X: Knots to be inserted
    :param r: len(x)
    :param Ubar: New Knots
    :param Qw: New Control Pts
    :return: None
    """
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

    for j in range(r, -1, -1):
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


def DegreeElevateCurve(n, p, U, Pw, t):
    """
    将曲线的次数升高t次
    :param n: NURBS曲线控制点数量=n+1
    :param p: NURBS曲线次数
    :param U: Original Knots
    :param Pw: Original Control Pts
    :param t:  level to be promoted
    :return: nh, Uh, Qw
    """

    nh = 0
    Uh = np.zeros(1000)
    Qw = np.zeros(1000)

    m = n + p + 1
    ph = p + t
    ph2 = ph / 2

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


def merge(lhs, rhs):
    i = 0
    j = 0
    ret = []
    while i != len(lhs) and j != len(rhs):
        if lhs[i] != rhs[j]:
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


def GetDistance(dim, lhs, rhs):
    assert len(lhs) == len(rhs) == dim

    l = 0
    for i in range(0, dim):
        l += math.pow(lhs[i] - rhs[i], 2)

    return math.sqrt(l)


class Skining(object):
    def __init__(self, _crvs):
        self.crv = _crvs
        self.common_knots = []

    def promote(self, p):
        p_list = []
        for crv in self.crv:
            p_list.append(crv.p)
        cpm = np.amax(p_list)

        if cpm > p:
            raise ValueError("Target order is less than existing max order!")

        crv_param_list = []
        for crv in self.crv:
            if p != crv.p:
                crv_param_list.append(DegreeElevateCurve(crv.n, crv.p, crv.knots, crv.ctrl_pts, p - crv.p))
            else:
                crv_param_list.append([crv.n, crv.knots, crv.ctrl_pts])

        return crv_param_list

    def merge_all(self):
        if len(self.crv) == 0:
            raise Exception("Empty Curve container!")

        ccu = self.crv[0].knots

        for i in range(1, len(self.crv)):
            ccu = merge(ccu, self.crv[i].knots)

        self.common_knots = ccu

    def calc_v_knots(self, q):
        K = len(self.crv)
        N = self.crv[0].n + 1
        d = np.zeros(N, float)
        seg_len = np.zeros((N, K - 1), float)

        for i in range(0, N):
            for j in range(1, K):
                seg_len[i][j - 1] = GetDistance(self.crv[i].ctrl_pts[j].shape[0], self.crv[i].ctrl_pts[j], self.crv[i].ctrl_pts[j - 1])
                d[i] += seg_len[i][j - 1]

        v = np.zeros(K, float)
        v[K - 1] = 1.0
        for k in range(1, K - 2):
            tmp = 0.0
            for i in range(0, N):
                tmp += seg_len[i][k - 1] / d[i]
            v[k] = v[k - 1] + tmp / N

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

        return v, knots

    def generate(self, p=5, q=5):
        self.promote(p)
        self.merge_all()
        self.calc_v_knots(q)


class Surface(object):
    def __init__(self):
        pass
