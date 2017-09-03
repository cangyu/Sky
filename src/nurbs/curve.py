from copy import deepcopy
from numpy.linalg import norm
from scipy.linalg import solve
from scipy.interpolate import BSpline, make_interp_spline
from scipy.integrate import romberg
from scipy.misc import comb
from src.nurbs.utility import *
from src.nurbs.transform import DCM, Quaternion
from src.iges.iges_entity110 import IGES_Entity110
from src.iges.iges_entity126 import IGES_Entity126


class ClampedNURBSCrv(object):
    def __init__(self, u, pw):
        """
        NURBS曲线
        :param u:节点矢量 
        :param pw: 带权控制点
        """

        self.U = np.copy(u)
        self.Pw = np.copy(pw)

        self.spl = BSpline(self.U, self.Pw, self.p)

    @property
    def m(self):
        """
        最后一个节点下标
        """

        return len(self.U) - 1

    @property
    def n(self):
        """
        最后一个控制点下标
        """

        return len(self.Pw) - 1

    @property
    def p(self):
        """
        曲线次数
        """

        return self.m - self.n - 1

    @property
    def weight(self):
        """
        权系数
        """

        return self.Pw[:, -1]

    @property
    def cpt(self):
        """
        求控制点
        """

        tn = len(self.Pw)
        ans = np.zeros((tn, 3))
        for i in range(tn):
            ans[i] = to_cartesian(self.Pw[i])
        return ans

    def __call__(self, u, d=0, return_cartesian=True):
        """
        计算曲线上的点
        :param u: 目标参数
        :param d: 求导次数
        :param return_cartesian: 返回结果是否带权
        :return: 曲线在u处的d阶导矢
        """

        pw = self.spl(u, d)
        return to_cartesian(pw) if return_cartesian else pw

    def length(self):
        """
        近似计算曲线长度
        """

        return romberg(lambda u: norm(self.__call__(u, 1)), self.U[0], self.U[-1])

    def curvature(self, u):
        """
        计算在给定位置处曲率
        """

        p1 = self.__call__(u, 1)
        p2 = self.__call__(u, 2)

        dd = np.zeros(3)
        dd[0] = math.pow(p2[2] * p1[1] - p2[1] * p1[2], 2)
        dd[1] = math.pow(p2[0] * p1[2] - p2[2] * p1[0], 2)
        dd[2] = math.pow(p2[1] * p1[0] - p2[0] * p1[1], 2)
        dividend = math.sqrt(np.sum(dd))

        dv = np.zeros(3)
        dv[0] = math.pow(p1[0], 2)
        dv[1] = math.pow(p1[1], 2)
        dv[2] = math.pow(p1[2], 2)
        divisor = math.pow(np.sum(dv), 3 / 2)

        kappa = dividend / divisor
        return kappa

    def reset(self, u, pw):
        """
        根据新的节点矢量和带权控制点重新初始化曲线
        :param u: New Knot vector.
        :param pw: New control points.
        :return: None.
        """

        self.U = np.copy(u)
        self.Pw = np.copy(pw)
        self.spl = BSpline(self.U, self.Pw, self.p)

    def reverse(self):
        """
        曲线反向
        """

        nu = np.full(self.m + 1, self.U[0] + self.U[-1]) - self.U[::-1]
        npw = self.Pw[::-1, :]
        self.reset(nu, npw)

    def pan(self, delta):
        """
        曲线整体平移
        :param delta: 偏移矢量
        :return: None.
        """

        dv = np.zeros(3)
        array_smart_copy(delta, dv)
        npw = np.empty(self.Pw.shape, float)

        for i in range(self.n + 1):
            cv = to_cartesian(self.Pw[i])
            cv += dv
            npw[i] = to_homogeneous(cv, self.Pw[i][-1])

        self.reset(self.U, npw)

    def rotate(self, ref, ax, ang):
        """
        将曲线绕过指定点的转轴旋转一定角度
        :param ref: 参考点
        :param ax: 旋转轴方向向量，按右手定则确定角度的正方向
        :param ang: 旋转角(Degree)
        :return: None
        """

        q = Quaternion.from_u_theta(ax, math.radians(ang))
        npw = np.empty(self.Pw.shape, float)
        for i in range(self.n + 1):
            cv = to_cartesian(self.Pw[i])
            cv -= ref
            cv = ref + q.rotate(cv)
            npw[i] = to_homogeneous(cv, self.Pw[i][-1])

        self.reset(self.U, npw)

    @property
    def start(self):
        """
        曲线起点
        """

        return to_cartesian(self.Pw[0])

    @property
    def end(self):
        """
        曲线终点
        """

        return to_cartesian(self.Pw[-1])

    def to_iges(self, *args, **kwargs):
        """
        将曲线以IGES标准中的第126号实体呈现
        :param args: 依次指示:isPlaner, isPeriodic, norm, form
        :param kwargs: Not used currently.
        :return: IGES_Entity126 Object.
        """

        planar = 0
        periodic = 0
        norm_vector = np.zeros(3)
        form = kwargs['form'] if 'form' in kwargs else 0
        if len(args) != 0:
            planar = args[0]
            periodic = args[1]
            norm_vector = args[2]

        w = self.weight
        cpt = self.cpt
        poly = 0 if (w != np.ones(w.shape)).any() else 1
        closed = 1 if equal(norm(self.end - self.start), 0.0) else 0

        return IGES_Entity126(self.p, self.n, planar, closed, poly, periodic, self.U, w, cpt, self.U[0], self.U[-1], norm_vector, form)

    def insert_knot(self, u, r=1):
        """
        插入一个节点若干次
        :param u: 待插入节点
        :param r: 插入的次数，要求s+r<=p, 其中s为u在原节点矢量中的重复度,p为曲线次数
        :return: None.
        """

        if r < 0:
            raise AssertionError('Invalid times!')

        if r == 0:
            return

        '''Insert'''
        s = sum(x == u for x in self.U)  # Counts of duplicates
        if s + r > self.p:
            raise ValueError('Too many Knot: {}\nExisting: {}, Targeting: {} Max: {}'.format(u, s, s + r, self.p))

        k = find_span(self.n, self.p, u, self.U)
        nU = np.insert(self.U, k + 1, np.full(r, u, float))  # New knot vector
        nPw = np.zeros((self.n + r + 1, 4))  # New homogeneous control points

        '''Calculate new control points'''
        Rw = np.zeros((self.p + 1, 4))  # Holding temporary points

        '''Store unchanged control points'''
        for i in range(k - self.p + 1):
            nPw[i] = np.copy(self.Pw[i])
        for i in range(k - s, self.n + 1):
            nPw[i + r] = np.copy(self.Pw[i])
        for i in range(self.p - s + 1):
            Rw[i] = np.copy(self.Pw[k - self.p + i])

        '''Insert target knot r times'''
        L = 0
        for j in range(1, r + 1):
            L = k - self.p + j
            for i in range(self.p - j - s + 1):
                alpha = (u - self.U[L + i]) / (self.U[i + k + 1] - self.U[L + i])
                Rw[i] = alpha * Rw[i + 1] + (1.0 - alpha) * Rw[i]
            nPw[L] = np.copy(Rw[0])
            nPw[k + r - j - s] = np.copy(Rw[self.p - j - s])

        '''Load remaining control points'''
        for i in range(L + 1, k - s):
            nPw[i] = np.copy(Rw[i - L])

        '''Update'''
        self.reset(nU, nPw)

    def refine(self, extra_knots):
        """
        节点细化，插入额外的节点序列
        :param extra_knots: 待插入节点序列(已按升序排好)
        """

        if len(extra_knots) == 0:
            return

        r = len(extra_knots) - 1
        nu = np.zeros(self.m + r + 2, float)  # New knot vector
        npw = np.zeros((self.n + r + 2, 4), float)  # New homogeneous control points

        '''Knot span'''
        a = find_span(self.n, self.p, extra_knots[0], self.U)
        b = find_span(self.n, self.p, extra_knots[r], self.U) + 1

        '''Copy unchanged control points and knots'''
        for j in range(a - self.p + 1):
            npw[j] = np.copy(self.Pw[j])
        for j in range(b - 1, self.n + 1):
            npw[j + r + 1] = np.copy(self.Pw[j])

        for j in range(a + 1):
            nu[j] = self.U[j]
        for j in range(b + self.p, self.m + 1):
            nu[j + r + 1] = self.U[j]

        '''Insert'''
        i = b + self.p - 1
        k = b + self.p + r
        for j in range(r, -1, -1):
            while extra_knots[j] <= self.U[i] and i > a:
                npw[k - self.p - 1] = np.copy(self.Pw[i - self.p - 1])
                nu[k] = self.U[i]
                k -= 1
                i -= 1

            npw[k - self.p - 1] = np.copy(npw[k - self.p])

            for l in range(1, self.p + 1):
                index = k - self.p + l
                alpha = nu[k + l] - extra_knots[j]
                if equal(alpha, 0.0):
                    npw[index - 1] = np.copy(npw[index])
                else:
                    alpha /= (nu[k + l] - self.U[i - self.p + l])
                    npw[index - 1] = alpha * npw[index - 1] + (1.0 - alpha) * npw[index]

            nu[k] = extra_knots[j]
            k -= 1

        self.reset(nu, npw)

    def remove_knot(self, u, num, delta=1e-6):
        """
        Remove specified knot 'num' times.
        :param u: Knot to be removed
        :type u: float
        :param num: Times of removal
        :type num: int
        :param delta: Max expected derivation
        :type delta: float
        :return: Times of actual removal
        :rtype: int
        """

        '''Defensive check'''
        if not (u in self.U):
            raise ValueError("Target knot not exist.")
        if equal(u, 0) or equal(u, 1):
            raise ValueError("Invalid input.")

        '''Find position and duplication'''
        r = 0
        while not equal(self.U[r], u):
            r += 1

        s = 0
        while equal(self.U[r], u):
            s += 1
            r += 1
        r -= 1

        '''Tolerance'''
        TOL = math.fabs(delta * min(self.weight) / (1 + max(list(map(lambda _p: norm(_p), self.cpt)))))

        '''Basic variables'''
        p = self.p
        m = self.m
        n = self.n
        ord = p + 1
        fout = (2 * r - s - p) // 2
        last = r - s
        first = r - p

        '''Temp'''
        temp = np.empty((2 * p + 1, 4), float)

        '''Removal'''
        t = 0
        while t < num:
            off = first - 1
            temp[0] = self.Pw[off]
            temp[last + 1 - off] = self.Pw[last + 1]
            i = first
            j = last
            ii = 1
            jj = last - off
            while j - i > t:
                alfi = (u - self.U[i]) / (self.U[i + ord + t] - self.U[i])
                alfj = (u - self.U[j - t]) / (self.U[j + ord] - self.U[j - t])
                temp[ii] = (self.Pw[i] - (1 - alfi) * temp[ii - 1]) / alfi
                temp[jj] = (self.Pw[j] - alfj * temp[jj + 1]) / (1 - alfj)
                i += 1
                ii += 1
                j -= 1
                jj -= 1

            if j - i < t:
                remflag = pnt_dist(temp[ii - 1], temp[jj + 1]) <= TOL
            else:
                alfi = (u - self.U[i]) / (self.U[i + ord + t] - self.U[i])
                tpnt = alfi * temp[ii + t + 1] + (1 - alfi) * temp[ii - 1]
                remflag = pnt_dist(self.Pw[i], tpnt) <= TOL

            if not remflag:
                break
            else:
                i = first
                j = last
                while j - i > t:
                    self.Pw[i] = temp[i - off]
                    self.Pw[j] = temp[j - off]
                    i += 1
                    j -= 1

            first -= 1
            last += 1
            t += 1

        if t == 0:
            return t

        for k in range(r + 1, m + 1):
            self.U[k - t] = self.U[k]

        j = fout
        i = j
        for k in range(1, t):
            if k % 2 == 0:
                j -= 1
            else:
                i += 1

        for k in range(i + 1, n + 1):
            self.Pw[j] = self.Pw[k]
            j += 1

        '''Drop tailing knot and control point'''
        if t != 0:
            self.reset(self.U[:-t], self.Pw[:-t])

        return t

    @classmethod
    def decompose(cls, crv):
        """
        Decompose the NURBS curve into several bezier segments.
        This is knot insertion in essence, just on its intrinsic knots.
        Optimization are performed especially.
        :param crv: Curve to be decomposed
        :type crv: ClampedNURBSCrv
        :return: Bezier segments
        """

        '''New knot vector and control points'''
        val = np.unique(crv.U)
        sorted(val)
        Qw = np.empty((len(val) - 1, crv.p + 1, 4), float)

        '''Calculate new control points'''
        alphas = np.empty(crv.p, float)
        a = crv.p
        b = crv.p + 1
        nb = 0
        for i in range(crv.p + 1):
            Qw[nb][i] = np.copy(crv.Pw[i])

        while b < crv.m:
            i = b
            while b < crv.m and equal(crv.U[b + 1], crv.U[b]):
                b += 1
            mult = b - i + 1
            if mult < crv.p:
                numer = crv.U[b] - crv.U[a]
                j = crv.p
                while j > mult:
                    alphas[j - mult - 1] = numer / (crv.U[a + j] - crv.U[a])
                    j -= 1
                r = crv.p - mult
                for j in range(1, r + 1):
                    save = r - j
                    s = mult + j
                    k = crv.p
                    while k >= s:
                        alpha = alphas[k - s]
                        Qw[nb][k] = alpha * Qw[nb][k] + (1.0 - alpha) * Qw[nb][k - 1]
                        k -= 1
                    if b < crv.m:
                        Qw[nb + 1][save] = np.copy(Qw[nb][crv.p])

            nb += 1
            if b < crv.m:
                for i in range(crv.p - mult, crv.p + 1):
                    Qw[nb][i] = np.copy(crv.Pw[b - crv.p + i])
                a = b
                b += 1

        '''Defensive Check'''
        if nb != len(Qw):
            raise AssertionError("Internal Error.")

        ret = []
        kidx = 0
        for i in range(nb):
            crv = BezierCrv(val[kidx], val[kidx + 1], crv.p, Qw[i])
            ret.append(crv)
            kidx += 1

        return ret

    def elevate(self, t):
        """
        将曲线升阶t次
        :param t: 升阶次数
        :type t: int
        :return: None
        """

        if t <= 0:
            return

        p = self.p
        val, cnt = np.unique(self.U, return_counts=True)

        '''Decompose'''
        bezier_seg = ClampedNURBSCrv.decompose(self)

        '''Merge with degree elevation'''
        new_crv = ClampedNURBSCrv.merge(bezier_seg, p + t)

        '''Elimination'''
        for k, u in enumerate(val):
            rts = p - cnt[k]
            if rts > 0:
                arts = new_crv.remove_knot(u, rts)
                if arts != rts:
                    raise RuntimeError("Failed to eliminate knot {}".format(u))

        '''Update'''
        self.reset(new_crv.U, new_crv.Pw)

    def reparameterization(self, alpha, beta, gamma, delta):
        """
        使用线性有理函数进行重新参数化
        不需要改变控制点，但权系数需要改变
        Denote 'u' as current parameter, 's' as target parameter,
        relations between 's' and 'u' are given as follows:
            alpha*u + beta
        s = ---------------
            gamma*u + delta
        :param alpha: Linear rational reparameterization parameter
        :type alpha: float
        :param beta: Linear rational reparameterization parameter
        :type beta: float
        :param gamma: Linear rational reparameterization parameter
        :type gamma: float
        :param delta: Linear rational reparameterization parameter
        :type delta: float
        :return: None
        """

        if not alpha * delta - gamma * beta > 0:
            raise AssertionError("Bad reparameterization.")

        def mu(u):
            return gamma * u + delta

        '''Calculate new knot'''
        old_knot = self.U
        new_knot = np.copy(list(map(lambda u: (alpha * u + beta) / (gamma * u + delta), old_knot)))

        '''Calculate new weight'''
        old_wt = self.weight
        new_wt = np.empty_like(old_wt)

        cp = self.p
        factor = 1.0
        for i in range(cp):
            factor *= mu(old_knot[i])

        wn = len(old_wt)
        for i in range(wn):
            factor /= mu(old_knot[i])
            factor *= mu(old_knot[i + cp])
            new_wt[i] = old_wt[i] / factor

        '''Calculate new weighted-control-pts'''
        cur_cpt = self.cpt
        wpt = np.empty((wn, 4), float)
        for i in range(wn):
            wpt[i] = to_homogeneous(cur_cpt[i], new_wt[i])

        '''Update'''
        self.reset(new_knot, wpt)

    def standard_reparameterization(self):
        """
        通过线性变换，将节点转为[0, 1]上的标准参数化
        :return: None
        """

        a = self.U[0]
        b = self.U[-1]
        alpha = 1.0
        beta = -a
        gamma = 0.0
        delta = b - a

        self.reparameterization(alpha, beta, gamma, delta)

    @classmethod
    def split(cls, crv, break_pts):
        """
        Split the curve into several segments.
        :param crv: NURBS curve to be split
        :type crv: ClampedNURBSCrv
        :param break_pts: Splitting knots
        :return: Curve segments
        """

        '''Count current knots'''
        val, cnt = np.unique(crv.U, return_counts=True)
        knot_dict = dict(zip(val, cnt))

        '''Calculate knots to be inserted'''
        sbp = sorted(break_pts)
        if sbp[0] <= crv.U[0] or sbp[-1] >= crv.U[-1]:
            raise ValueError("Invalid break knot.")

        bkt = []
        cp = crv.p
        for u in sbp:
            exist_cnt = knot_dict.get(u) if u in knot_dict else 0
            tc = cp - exist_cnt
            for k in range(tc):
                bkt.append(u)

        '''Insert breaking knots'''
        crv0 = deepcopy(crv)
        crv0.refine(np.copy(bkt))

        '''Extract each segment'''
        ret = []
        sbp.append(crv0.U[-1])
        knot_num = crv0.m + 1
        prev_knot = crv0.U[0]
        cki = cp + 1
        cpi = 0
        for u in sbp:
            ck = []
            for i in range(cp + 1):
                ck.append(prev_knot)
            while cki < knot_num and crv0.U[cki] <= u:
                ck.append(crv0.U[cki])
                cki += 1
            if cki < knot_num:
                ck.append(u)
            prev_knot = u
            cpn = len(ck) - cp - 1
            cpi_next = cpi + cpn
            csg = ClampedNURBSCrv(ck, crv0.Pw[cpi:cpi_next])
            csg.standard_reparameterization()
            ret.append(csg)
            cpi = cpi_next - 1

        return ret

    @classmethod
    def merge(cls, crv_list, p=None):
        """
        Merge several bezier curves into one NURBS curve.
        :param crv_list: Bezier curve list
        :param p: Target order
        :type p: int
        :return: Merged curve
        :rtype: ClampedNURBSCrv
        """

        '''Do not affect original data'''
        bezier_list = deepcopy(crv_list)

        '''Check continuity'''
        prev_ending = bezier_list[0].start
        for crv in bezier_list:
            if not isinstance(crv, BezierCrv):
                raise AssertionError("Invalid input.")
            if not equal(norm(crv.start - prev_ending), 0):
                raise AssertionError("Not continuous.")
            prev_ending = crv.end

        '''Check Order'''
        crv_order = 0 if p is None else p
        for crv in bezier_list:
            crv_order = max(crv_order, crv.p)

        '''Degree elevation'''
        for k, crv in enumerate(bezier_list):
            t = crv_order - crv.p
            if t > 0:
                bezier_list[k].elevate(t)

        '''Assembly'''
        return bezier_merge(bezier_list)


class BezierCrv(ClampedNURBSCrv):
    def __init__(self, a, b, p, pw):
        kv = []
        for i in range(p + 1):
            kv.append(a)
        for i in range(p + 1):
            kv.append(b)

        super(BezierCrv, self).__init__(kv, pw)

    @property
    def a(self):
        return self.U[0]

    @property
    def b(self):
        return self.U[-1]

    def elevate(self, t):
        """
        将Bezier曲线升阶t次
        :param t: 升阶次数
        :type t: int
        :return: None.
        """

        bezier_deg_elev(self, t)


def bezier_deg_elev(crv, t):
    """
    Degree elevation of an Bezier curve
    :param crv: Bezier curve to be elevated
    :type crv: BezierCrv
    :param t: Elevation level
    :type t: int
    :return: None
    """

    if t <= 0:
        return

    nh = ph = crv.p + t + 1
    npw = np.zeros((nh, 4))
    kv = []

    '''Knots'''
    for i in range(ph):
        kv.append(crv.a)
    for i in range(ph):
        kv.append(crv.b)

    '''Control points'''
    for i in range(nh):
        for j in range(max(0, i - t), min(crv.p, i) + 1):
            cc = comb(crv.p, j, exact=True) * comb(t, i - j, exact=True) / comb(crv.p + t, i, exact=True)
            npw[i] += cc * crv.Pw[j]

    '''Update'''
    crv.reset(kv, npw)


def bezier_merge(bezier_list):
    """
    Merge a set of bezier curves.
    We assume that the input curve set is continous and have same order.
    :param bezier_list: A set of bezier curves to be merge
    :return: Curve with eliminated knots
    :rtype: ClampedNURBSCrv
    """

    crv_order = bezier_list[0].p
    seg_num = len(bezier_list)

    '''Construct knots'''
    nU = np.empty((seg_num + 1) * crv_order + 2, float)
    nU[0] = bezier_list[0].a
    k = 1
    for bsg in bezier_list:
        tmp = bsg.a
        for i in range(crv_order):
            nU[k] = tmp
            k += 1
    tmp = bezier_list[-1].b
    for i in range(crv_order + 1):
        nU[k] = tmp
        k += 1

    '''Construct control points'''
    nPw = np.empty((seg_num * crv_order + 1, 4), float)
    k = 0
    for bsg in bezier_list:
        for i in range(bsg.n):
            nPw[k] = np.copy(bsg.Pw[i])
            k += 1
    nPw[-1] = np.copy(bezier_list[-1].Pw[-1])

    '''Construct NURBS curve'''
    return ClampedNURBSCrv(nU, nPw)


class GlobalInterpolatedCrv(ClampedNURBSCrv):
    def __init__(self, pts, p=3, method='centripetal'):
        """
        构造一条p次非有理B样条曲线插值于pts
        :param pts: 待插值点序列
        :param p: 目标曲线次数
        :param method: 计算插值点参数的方法
        """

        n, dim = pts.shape
        n -= 1
        param = calc_pnt_param(pts, method)
        kv = calc_knot_vector(param, p)
        cpt = calc_ctrl_pts(kv, p, pts, param)

        pw = np.zeros((n + 1, dim + 1))
        for i in range(0, n + 1):
            pw[i] = to_homogeneous(cpt[i])

        super(GlobalInterpolatedCrv, self).__init__(kv, pw)


def calc_pnt_param(pts, method):
    """
    计算每个插值点所对应的参数。
    :param pts: 插值点坐标序列
    :param method: 参数计算方法
    :return: 插值点坐标序列对应的参数序列([0,1])
    """

    if method not in ['chord', 'centripetal']:
        raise ValueError("Invalid method parameter!")

    n = len(pts) - 1
    param = np.zeros(n + 1)
    param[n] = 1.0

    dist = np.zeros(n + 1)
    for i in range(1, n + 1):
        dist[i] = pnt_dist(pts[i - 1], pts[i])

    d = 0
    if method == 'chord':  # 弦长参数化
        for i in range(1, n + 1):
            d += dist[i]
    else:  # 向心参数化，数据点急转弯变化时效果好
        for i in range(1, n + 1):
            dist[i] = math.sqrt(dist[i])
            d += dist[i]

    for i in range(1, n):
        param[i] = param[i - 1] + dist[i] / d

    return param


def calc_knot_vector(param, p):
    """
    取平均值方法计算节点
    :param param: 插值点序列对应的参数序列
    :param p: 目标曲线次数
    :return: 目标曲线节点矢量([0,1])
    """

    n = len(param) - 1
    m = n + p + 1
    knots = np.zeros(m + 1)

    '''Tail'''
    for i in range(0, p + 1):
        knots[m - i] = 1.0

    '''Prepare'''
    acc = 0.0
    for i in range(0, p):
        acc += param[i]

    '''Iterate'''
    for j in range(1, n - p + 1):
        acc -= param[j - 1]
        acc += param[p - 1 + j]
        knots[p + j] = acc / p

    return knots


def calc_ctrl_pts(U, p, pts, param):
    """
    求解线性方程组得到控制点
    :param U: 节点矢量
    :param p: 目标曲线次数
    :param pts: 插值点序列
    :param param: 插值点所对应参数
    :return: 控制点序列
    """
    n, dim = pts.shape
    n -= 1

    ctrl_pts = np.zeros((n + 1, dim))

    '''Coefficient Matrix'''
    cm = np.zeros((n + 1, n + 1))
    for k in range(0, n + 1):
        cm[k] = all_basis_val(param[k], p, U)

    '''Solve'''
    Q = np.zeros((dim, n + 1))
    P = np.zeros((dim, n + 1))

    for i in range(0, dim):
        for j in range(0, n + 1):
            Q[i][j] = pts[j][i]

    for i in range(0, dim):
        P[i] = solve(cm, Q[i])

    for i in range(0, n + 1):
        for j in range(0, dim):
            ctrl_pts[i][j] = P[j][i]

    return ctrl_pts


class Spline(ClampedNURBSCrv):
    def __init__(self, pts, p=3, bc=([(2, (0, 0, 0))], [(2, (0, 0, 0))]), method='centripetal'):
        """
        带端点切矢量的全局曲线插值
        Note:
        此处也可以将bc取为None，从而和GlobalInterpolatedCrv功能相同,
        但SciPy中的默认参数化方法可能不一样, 经测试其构造knot的方法可能不是简单地取平均，有待考证
        :param pts: 待插值数据点
        :param p: 插值曲线次数
        :type p: int
        :param bc: 在两端点处的边界条件，默认取自然边界条件
        """

        sp = calc_pnt_param(np.copy(pts), method)
        f = make_interp_spline(sp, np.copy(pts), k=p, bc_type=bc)
        pw = np.ones((len(f.t) - p - 1, 4), float)
        for k, pnt in enumerate(f.c):
            for d in range(3):
                pw[k][d] = pnt[d]

        super(Spline, self).__init__(f.t, pw)


class Line(ClampedNURBSCrv):
    def __init__(self, a, b):
        """
        两点间直线段
        :param a: 起始点坐标
        :param b: 终点坐标
        """

        u = np.array([0, 0, 1, 1])
        pw = np.array([[0, 0, 0, 1], [0, 0, 0, 1]], float)
        array_smart_copy(a, pw[0])
        array_smart_copy(b, pw[1])

        super(Line, self).__init__(u, pw)

    def length(self):
        return pnt_dist(self.start, self.end)

    def curvature(self, u):
        return 0.0

    def to_iges(self):
        return IGES_Entity110(to_cartesian(self.Pw[0]), to_cartesian(self.Pw[-1]))


class Arc(ClampedNURBSCrv):
    def __init__(self, r, theta):
        """
        XY平面内简化圆弧，以原点为圆心，起始点为(r,0),法向量为(0,0,1)
        :param r: 半径
        :param theta: 圆心角, will be fitted into range (0,360]
        """

        while theta <= 0:
            theta += 360
        while theta > 360:
            theta -= 360

        self.radius = r
        self.theta = theta

        narcs = int(np.ceil(theta / 90))
        theta = np.deg2rad(theta)
        dtheta = theta / narcs
        w1 = np.cos(dtheta / 2)
        dknot = 1.0 / narcs

        m = 2 * narcs + 3  # 最后一个节点下标
        n = 2 * narcs  # 最后一个控制点下标

        U = np.zeros(m + 1)
        P = np.zeros((n + 1, 3))
        Pw = np.zeros((n + 1, 4))

        '''Knot Vector'''
        U[-1] = U[-2] = U[-3] = 1.0
        for i in range(1, narcs):
            cur_index = 1 + 2 * i
            U[cur_index] = U[cur_index + 1] = i * dknot

        '''Control Points'''
        P[0] = np.array([r, 0, 0], float)
        Pw[0] = to_homogeneous(P[0], 1.0)
        T0 = np.array([0.0, 1.0, 0.0])

        index = 0
        angle = 0.0
        for i in range(1, narcs + 1):
            angle += dtheta
            P[index + 2] = np.array([r * np.cos(angle), r * np.sin(angle), 0.0])
            Pw[index + 2] = to_homogeneous(P[index + 2], 1.0)
            T2 = np.array([-np.sin(angle), np.cos(angle), 0.0])
            P[index + 1] = line_intersection(P[index], T0, P[index + 2], T2)
            Pw[index + 1] = to_homogeneous(P[index + 1], w1)
            index += 2
            if i < narcs:
                T0 = T2

        super(Arc, self).__init__(U, Pw)

    def curvature(self, u):
        return 1.0 / self.radius

    def length(self):
        return self.radius * np.deg2rad(self.theta)

    @classmethod
    def from_2pnt(cls, start_pnt, end_pnt, theta, norm_vector):
        """
        空间圆弧
        :param start_pnt: 起始点坐标
        :param end_pnt: 终止点坐标
        :param theta: 圆心角
        :param norm_vector: 所在平面的法向量，按右手定则， 四指依次扫过start_pnt和end_pnt
        """

        '''Basic Variables'''
        sp = np.copy(start_pnt)
        ep = np.copy(end_pnt)
        theta = np.deg2rad(theta)
        radius = 0.5 * pnt_dist(sp, ep) / np.sin(theta / 2)
        w = radius * np.cos(theta / 2)
        cdir = normalize(np.cross(norm_vector, ep - sp))
        center = 0.5 * (sp + ep) + cdir * w

        '''Rotate and pan'''
        arc = cls(radius, np.rad2deg(theta))
        base1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
        xdir = sp - center
        base2 = np.array([xdir, np.cross(norm_vector, xdir), norm_vector], float)

        mrot = DCM(base1, base2).rot_matrix
        mrot_trans = np.transpose(mrot)
        pts = np.copy(arc.cpt)
        wg = np.copy(arc.weight)
        pw = np.zeros((len(pts), 4))
        for i in range(len(pts)):
            pts[i] = pts[i] * mrot_trans
            pts[i] += center
            pw[i] = to_homogeneous(pts[i], wg[i])

        '''Reconstruct'''
        arc.reset(arc.U, pw)
        return arc


class ConicArc(ClampedNURBSCrv):
    def __init__(self, _p0, _t0, _p2, _t2, _p):
        """
        单段有理Bezier圆锥截线弧
        :param _p0: 起始点
        :param _t0: 起始点处切矢量
        :param _p2: 终止点
        :param _t2: 终止点处切矢量
        :param _p: 曲线上一点坐标
        """

        p0 = np.copy(_p0)
        t0 = np.copy(_t0)
        p2 = np.copy(_p2)
        t2 = np.copy(_t2)
        p = np.copy(_p)

        '''Knots'''
        nu = np.array([0, 0, 0, 1, 1, 1], float)

        '''Calculate mid-pnt weight and coordinate'''
        v02 = p2 - p0
        if not np.cross(t0, t2).any():
            w1 = 0.0
            alf0, alf2, p1 = line_intersection(p, t0, p0, v02, True)
            a = math.sqrt(alf2 / (1 - alf2))
            u = a / (1 + a)
            b = 2 * u * (1 - u)
            b = -alf0 * (1 - b) / b
            p1 = b * t0
        else:
            p1 = line_intersection(p0, t0, p2, t2)
            v1p = p - p1
            alf0, alf2, tmp = line_intersection(p1, v1p, p0, v02, True)
            a = math.sqrt(alf2 / (1 - alf2))
            u = a / (1 + a)
            num = math.pow(1 - u, 2) * np.dot(p - p0, p1 - p) + math.pow(u, 2) * np.dot(p - p2, p1 - p)
            den = 2 * u * (1 - u) * np.dot(p1 - p, p1 - p)
            w1 = num / den

        npw = np.empty((3, 4), float)
        npw[0] = to_homogeneous(p0, 1)
        npw[1] = to_homogeneous(p1, w1)
        npw[2] = to_homogeneous(p2, 1)

        '''Set-up'''
        super(ConicArc, self).__init__(nu, npw)


class LocalCubicInterpolatedCrv(ClampedNURBSCrv):
    def __init__(self, pts, tv):
        """
        Interpolate points using cubic bezier curve
        with specified tangent vector on each segment ending.
        :param pts: Points to be interpolated.
        :param tv: Tangent vector.
        """

        if len(pts) != len(tv):
            raise AssertionError("Inconsistent input.")

        n = len(pts) - 1
        pw = np.ones((2 * n + 2, 4))
        u = np.zeros(2 * n + 6)

        '''Init'''
        array_smart_copy(pts[0], pw[0])
        array_smart_copy(pts[-1], pw[-1])
        u[-1] = u[-2] = 1.0
        utv = np.empty_like(tv)
        for k, t in enumerate(tv):
            utv[k] = normalize(t)

        '''Build Pw'''
        for i in range(n):
            t1 = utv[i + 1] + utv[i]
            t2 = pts[i + 1] - pts[i]
            a = 16 - norm(t1) ** 2
            b = 12 * np.dot(t1, t2)
            c = -36 * norm(t2) ** 2
            dlt = b ** 2 - 4 * a * c

            if dlt < 0:
                raise ValueError("Invalid delta.")

            dlt = math.sqrt(dlt)
            x0 = 0.5 * (-b + dlt) / a
            x1 = 0.5 * (-b - dlt) / a
            alpha = x0 if x0 >= 0 else x1

            p1 = pts[i] + alpha / 3 * utv[i]
            p2 = pts[i + 1] - alpha / 3 * utv[i + 1]
            array_smart_copy(p1, pw[2 * i + 1])
            array_smart_copy(p2, pw[2 * i + 2])

        '''Build Knot'''
        prev = 0
        for i in range(1, n + 1):
            k = i - 1
            cur = prev + 3 * pnt_dist(pw[2 * k + 1][:3], pts[k])
            prev = cur
            u[2 + 2 * i] = u[3 + 2 * i] = cur

        for i in range(4, len(u) - 2):
            u[i] /= prev

        super(LocalCubicInterpolatedCrv, self).__init__(u, pw)
