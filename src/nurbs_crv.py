import unittest
import math
import sys
from copy import deepcopy
import numpy as np
from numpy.linalg import norm
from scipy.integrate import romberg
from scipy.interpolate import BSpline, make_interp_spline
from scipy.linalg import solve
from scipy.misc import comb
from iges import Model, Entity110, Entity126
from transform import Quaternion, DCM
from misc import array_smart_copy, normalize, pnt_dist, read_airfoil_pts, sqrt2, sqrt3
from nurbs_basis import to_homogeneous, to_cartesian, all_basis_val, find_span, line_intersection

"""
Implementation of the NURBS curve.

Note:
All the NURBS notations are in the 'Clamped' format by default.

TODO:
(1) Optimize the calculation of derivatives inside 'Crv.__call__(u, d)'.
(2) Determine better converge criteria inside 'point_inverse(c, p)'.
"""


class Crv(object):
    def __init__(self, u, pw):
        """
        NURBS Curve
        :param u: Knot vector.
        :param pw: Control points with weights(in homogeneous format).
        """

        self.U = np.copy(u)
        self.Pw = np.copy(pw)

        self.spl = BSpline(self.U, self.Pw, self.p)

    def __repr__(self):
        ret = 'NURBS Curve in clamped format.\n'
        ret += 'Knot Vector:\n{}\n'.format(self.U)
        ret += 'Ctrl point:\n{}\n'.format(self.cpt)
        ret += 'Weight:\n{}\n'.format(self.weight)
        return ret

    @property
    def m(self):
        """
        The last index of knots.
        """

        return len(self.U) - 1

    @property
    def n(self):
        """
        The last index of control points.
        """

        return len(self.Pw) - 1

    @property
    def p(self):
        """
        Degree of the curve.
        """

        return self.m - self.n - 1

    @property
    def start(self):
        """
        Start of the curve.
        :return: Point in cartesian format.
        """

        return to_cartesian(self.Pw[0])

    @property
    def end(self):
        """
        End of the curve.
        :return: Point in cartesian format.
        """

        return to_cartesian(self.Pw[-1])

    @property
    def weight(self):
        """
        Get the weight sequence for all control points.
        """

        return self.Pw[:, -1]

    @property
    def cpt(self):
        """
        Get all the control points.
        """

        tn = len(self.Pw)
        ans = np.empty((tn, 3))
        for i in range(tn):
            ans[i] = to_cartesian(self.Pw[i])
        return ans

    def to_iges(self, *args, **kwargs):
        """
        Represent the curve in IGES_Entity126 format.
        :param args: Denote 'isPlaner', 'isPeriodic', 'norm' in sequence.
        :param kwargs: 'form' maybe denoted.
        :return: Curve in IGES_Entity126 format.
        :rtype: Entity126
        """

        form = kwargs['form'] if 'form' in kwargs else 0
        if len(args) != 0:
            planar = args[0]
            periodic = args[1]
            norm_vector = args[2]
        else:
            planar = 0
            periodic = 0
            norm_vector = np.zeros(3)

        w = self.weight
        cpt = self.cpt
        poly = 0 if (w != np.ones(w.shape)).any() else 1
        closed = 1 if math.isclose(norm(self.end - self.start), 0) else 0

        return Entity126(self.p, self.n, planar, closed, poly, periodic, self.U, w, cpt, self.U[0], self.U[-1], norm_vector, form)

    def __call__(self, u, d=0):
        """
        Calculate the point corresponding to given parameter.
        :param u: Target parameter.
        :param d: Degree of derivation.
        :type d: int
        :return: Value at u with d times derivation.
        """

        aw = np.copy(list(map(lambda der: self.spl(u, der), range(d + 1))))
        ck = np.empty((d + 1, 3), float)

        for k in range(d + 1):
            v = aw[k][:3]
            for i in range(1, k + 1):
                v -= comb(k, i) * aw[i][-1] * ck[k - i]
            ck[k] = v / aw[0][-1]

        return ck[d]

    @property
    def length(self):
        """
        Calculate the length of the whole curve approximately.
        """

        return romberg(lambda u: norm(self.__call__(u, 1)), self.U[0], self.U[-1])

    def curvature(self, u):
        """
        Calculate the curvature at given position.
        :param u: Target parameter.
        :return: Curvature at given position.
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
        Reset the curve with new knot vector and weighted control points.
        :param u: New Knot vector.
        :param pw: New control points.
        :return: None.
        """

        self.U = np.copy(u)
        self.Pw = np.copy(pw)
        self.spl = BSpline(self.U, self.Pw, self.p)

    def reverse(self):
        """
        Reverse the curve parameterization without changing the geometry.
        :return: None.
        """

        nu = np.full(self.m + 1, self.U[0] + self.U[-1]) - self.U[::-1]
        npw = self.Pw[::-1, :]
        self.reset(nu, npw)

    def pan(self, delta):
        """
        Pan the whole curve in specific offset.
        :param delta: The delta vector.
        :return: None.
        """

        dv = np.zeros(3)
        array_smart_copy(delta, dv)
        npw = list(map(lambda _p: to_homogeneous(to_cartesian(_p) + dv, _p[-1]), self.Pw))
        self.reset(self.U, npw)

    def rotate(self, ref, ax, ang):
        """
        Rotate the curve with specific angle along specific rotation axis.
        :param ref: Anchor point of the rotation axis.
        :param ax: Direction vector of the rotation axis(positive direction is given by the right-hand rule).
        :param ang: Rotation angle(in degree).
        :return: None
        """

        q = Quaternion.from_u_theta(ax, math.radians(ang))
        npw = list(map(lambda pnt: to_homogeneous(ref + q.rotate(to_cartesian(pnt) - ref), pnt[-1]), self.Pw))
        self.reset(self.U, npw)

    def insert_knot(self, u, r=1):
        """
        Insert a knot several times.
        :param u: Knot to be inserted.
        :type u: float
        :param r: Times of insertion.
                  It's required that s+r<=p, where s is the multiplicity inside the original knot vector,
                  and p is the degree of the curve.
        :type r: int
        :return: None.
        """

        if r < 0:
            raise ValueError('Invalid times!')
        if r == 0:
            return

        '''Insert'''
        s = sum(x == u for x in self.U)  # Counts of duplicates
        if s + r > self.p:
            raise ValueError('Too many Knot: {}, existing: {}, targeting: {}, max: {}.'.format(u, s, s + r, self.p))

        k = find_span(self.n, self.p, u, self.U)
        nu = np.insert(self.U, k + 1, np.full(r, u, float))  # New knot vector
        npw = np.zeros((self.n + r + 1, 4))  # New homogeneous control points

        '''Calculate new control points'''
        rw = np.zeros((self.p + 1, 4))  # Holding temporary points

        '''Store unchanged control points'''
        for i in range(k - self.p + 1):
            npw[i] = np.copy(self.Pw[i])
        for i in range(k - s, self.n + 1):
            npw[i + r] = np.copy(self.Pw[i])
        for i in range(self.p - s + 1):
            rw[i] = np.copy(self.Pw[k - self.p + i])

        '''Insert target knot r times'''
        ll = 0
        for j in range(1, r + 1):
            ll = k - self.p + j
            for i in range(self.p - j - s + 1):
                alpha = (u - self.U[ll + i]) / (self.U[i + k + 1] - self.U[ll + i])
                rw[i] = alpha * rw[i + 1] + (1.0 - alpha) * rw[i]
            npw[ll] = np.copy(rw[0])
            npw[k + r - j - s] = np.copy(rw[self.p - j - s])

        '''Load remaining control points'''
        for i in range(ll + 1, k - s):
            npw[i] = np.copy(rw[i - ll])

        '''Update'''
        self.reset(nu, npw)

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
                if math.isclose(alpha, 0.0):
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
        if math.isclose(u, 0) or math.isclose(u, 1):
            raise ValueError("Invalid input.")

        '''Find position and duplication'''
        r = 0
        while not math.isclose(self.U[r], u):
            r += 1

        s = 0
        while math.isclose(self.U[r], u):
            s += 1
            r += 1
        r -= 1

        '''Tolerance'''
        tolerance = math.fabs(delta * min(self.weight) / (1 + max(list(map(lambda _p: norm(_p), self.cpt)))))

        '''Basic variables'''
        p = self.p
        m = self.m
        n = self.n
        order = p + 1
        f_out = (2 * r - s - p) // 2
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
                alfi = (u - self.U[i]) / (self.U[i + order + t] - self.U[i])
                alfj = (u - self.U[j - t]) / (self.U[j + order] - self.U[j - t])
                temp[ii] = (self.Pw[i] - (1 - alfi) * temp[ii - 1]) / alfi
                temp[jj] = (self.Pw[j] - alfj * temp[jj + 1]) / (1 - alfj)
                i += 1
                ii += 1
                j -= 1
                jj -= 1

            if j - i < t:
                remflag = pnt_dist(temp[ii - 1], temp[jj + 1]) <= tolerance
            else:
                alfi = (u - self.U[i]) / (self.U[i + order + t] - self.U[i])
                tpnt = alfi * temp[ii + t + 1] + (1 - alfi) * temp[ii - 1]
                remflag = pnt_dist(self.Pw[i], tpnt) <= tolerance

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

        j = f_out
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
        :type crv: Crv
        :return: Bezier segments
        """

        '''New knot vector and control points'''
        val = np.unique(crv.U)
        sorted(val)
        qw = np.empty((len(val) - 1, crv.p + 1, 4), float)

        '''Calculate new control points'''
        alphas = np.empty(crv.p, float)
        a = crv.p
        b = crv.p + 1
        nb = 0
        for i in range(crv.p + 1):
            qw[nb][i] = np.copy(crv.Pw[i])

        while b < crv.m:
            i = b
            while b < crv.m and math.isclose(crv.U[b + 1], crv.U[b]):
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
                        qw[nb][k] = alpha * qw[nb][k] + (1.0 - alpha) * qw[nb][k - 1]
                        k -= 1
                    if b < crv.m:
                        qw[nb + 1][save] = np.copy(qw[nb][crv.p])

            nb += 1
            if b < crv.m:
                for i in range(crv.p - mult, crv.p + 1):
                    qw[nb][i] = np.copy(crv.Pw[b - crv.p + i])
                a = b
                b += 1

        '''Defensive Check'''
        if nb != len(qw):
            raise AssertionError("Internal Error.")

        ret = []
        kidx = 0
        for i in range(nb):
            crv = BezierCrv(val[kidx], val[kidx + 1], crv.p, qw[i])
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
        bezier_seg = Crv.decompose(self)

        '''Merge with degree elevation'''
        new_crv = Crv.merge(bezier_seg, p + t)

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
        :type crv: Crv
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
            csg = Crv(ck, crv0.Pw[cpi:cpi_next])
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
        :rtype: Crv
        """

        '''Do not affect original data'''
        bezier_list = deepcopy(crv_list)

        '''Check continuity'''
        prev_ending = bezier_list[0].start
        for crv in bezier_list:
            if not isinstance(crv, BezierCrv):
                raise AssertionError("Invalid input.")
            if not math.isclose(norm(crv.start - prev_ending), 0):
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


class BezierCrv(Crv):
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
    We assume that the input curve set is continous and share common degree.
    :param bezier_list: A set of bezier curves to be merge
    :return: Curve with eliminated knots
    :rtype: Crv
    """

    crv_order = bezier_list[0].p
    seg_num = len(bezier_list)

    '''Construct knots'''
    nu = np.empty((seg_num + 1) * crv_order + 2, float)
    nu[0] = bezier_list[0].a
    k = 1
    for bsg in bezier_list:
        tmp = bsg.a
        for i in range(crv_order):
            nu[k] = tmp
            k += 1
    tmp = bezier_list[-1].b
    for i in range(crv_order + 1):
        nu[k] = tmp
        k += 1

    '''Construct control points'''
    npw = np.empty((seg_num * crv_order + 1, 4), float)
    k = 0
    for bsg in bezier_list:
        for i in range(bsg.n):
            npw[k] = np.copy(bsg.Pw[i])
            k += 1
    npw[-1] = np.copy(bezier_list[-1].Pw[-1])

    '''Construct NURBS curve'''
    return Crv(nu, npw)


class GlobalInterpolatedCrv(Crv):
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
        for i in range(n + 1):
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
    for i in range(p + 1):
        knots[m - i] = 1.0

    '''Prepare'''
    acc = 0.0
    for i in range(p):
        acc += param[i]

    '''Iterate'''
    for j in range(1, n - p + 1):
        acc -= param[j - 1]
        acc += param[p - 1 + j]
        knots[p + j] = acc / p

    return knots


def calc_ctrl_pts(u_vec, p, pts, param):
    """
    求解线性方程组得到控制点
    :param u_vec: 节点矢量
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
    for k in range(n + 1):
        cm[k] = all_basis_val(param[k], p, u_vec)

    '''Solve'''
    bq = np.zeros((dim, n + 1))
    bp = np.zeros((dim, n + 1))

    for i in range(dim):
        for j in range(0, n + 1):
            bq[i][j] = pts[j][i]

    for i in range(dim):
        bp[i] = solve(cm, bq[i])

    for i in range(n + 1):
        for j in range(dim):
            ctrl_pts[i][j] = bp[j][i]

    return ctrl_pts


class Spline(Crv):
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


class Line(Crv):
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

    @property
    def length(self):
        return pnt_dist(self.start, self.end)

    def curvature(self, u):
        return 0.0

    def to_iges(self, *args, **kwargs):
        return Entity110(to_cartesian(self.Pw[0]), to_cartesian(self.Pw[-1]))


class Arc(Crv):
    def __init__(self, center, start, theta, norm_vec):
        """
        Spatial Arc.
        The direction of the arc is defined by the right-hand rule.
        :param center: Coordinate of the center.
        :param start: Coordinate of the starting point.
        :param theta: Central angle(in degree).
        :type theta: float
        :param norm_vec: The norm vector of the plane on which the arc is located.
        """

        self.center = np.copy(center)
        self.radius = pnt_dist(center, start)
        self.theta = theta
        self.norm_vec = np.copy(norm_vec)

        '''Basic Arc'''
        nu, npw = Arc.construct_planar(self.radius, self.theta)
        ncp = np.copy(list(map(lambda u: to_cartesian(u), npw)))
        n = len(ncp) - 1

        '''Pan and Rotate'''
        sp = np.copy(start)
        nx = sp - self.center
        base1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        base2 = np.array([nx, np.cross(self.norm_vec, nx), self.norm_vec])
        rot_mtx = np.transpose(DCM(base1, base2).rot_matrix)
        for i in range(n + 1):
            ncp[i] = ncp[i] * rot_mtx + self.center
            npw[i] = to_homogeneous(ncp[i], npw[i][-1])

        super(Arc, self).__init__(nu, npw)

    @classmethod
    def construct_planar(cls, radius, theta):
        """
        Arc on XY plane, with center at (0, 0, 0)
        :param radius: Radius of the arc.
        :type radius: float
        :param theta: Central angle of the arc.
        :type theta: float
        :return: Knot vector and weighted control points.
        """

        if radius <= 0:
            raise ValueError('Invalid radius.')
        if theta <= 0 or theta > 360:
            raise ValueError('Invalid central angle.')

        arc_num = int(math.ceil(theta / 90))
        theta_rad = math.radians(theta)
        dlt_theta = theta_rad / arc_num
        w1 = math.cos(dlt_theta / 2)
        dlt_knot = 1.0 / arc_num

        m = 2 * arc_num + 3  # Subscript of the last knot vector
        n = 2 * arc_num  # Subscript of the last control points

        nu = np.zeros(m + 1)
        ncp = np.zeros((n + 1, 3))
        npw = np.zeros((n + 1, 4))

        '''Knot Vector'''
        nu[-1] = nu[-2] = nu[-3] = 1.0
        for i in range(1, arc_num):
            cur_index = 1 + 2 * i
            nu[cur_index] = nu[cur_index + 1] = i * dlt_knot

        '''Control Points'''
        ncp[0] = np.array([radius, 0, 0], float)
        npw[0] = to_homogeneous(ncp[0], 1.0)
        t0 = np.array([0.0, 1.0, 0.0])

        idx = 0
        angle = 0.0
        for i in range(1, arc_num + 1):
            angle += dlt_theta
            ncp[idx + 2] = np.array([radius * math.cos(angle), radius * math.sin(angle), 0.0])
            npw[idx + 2] = to_homogeneous(ncp[idx + 2], 1.0)
            t2 = np.array([-math.sin(angle), math.cos(angle), 0.0])
            ncp[idx + 1] = line_intersection(ncp[idx], t0, ncp[idx + 2], t2)
            npw[idx + 1] = to_homogeneous(ncp[idx + 1], w1)
            idx += 2
            if i < arc_num:
                t0 = t2

        return nu, npw

    @classmethod
    def from_2pnt(cls, start, end, theta, norm_vec):
        sp = np.copy(start)
        ep = np.copy(end)
        theta_rad = math.radians(theta)
        radius = 0.5 * pnt_dist(sp, ep) / math.sin(theta_rad / 2)
        w = radius * math.cos(theta_rad / 2)
        center_axis = normalize(np.cross(norm_vec, ep - sp))
        center = 0.5 * (sp + ep) + center_axis * w
        return cls(center, sp, theta, norm_vec)

    @classmethod
    def closed_circle(cls, radius, center=(0, 0, 0), norm_vec=(0, 0, 1)):
        ref = np.copy(center)
        base1 = np.array([0, 0, 1])
        base2 = np.copy(norm_vec)
        rot_axis = np.cross(base1, base2)
        rot_angle = math.degrees(math.acos(np.dot(base1, base2) / (norm(base1) * norm(base2))))
        nu, npw = Arc.construct_planar(radius, 360)
        if not math.isclose(rot_angle, 0):
            q = Quaternion.from_u_theta(rot_axis, rot_angle)
            npw = np.copy(list(map(lambda p: to_homogeneous(q.rotate(to_cartesian(p)), p[-1]), npw)))
        npw = np.copy(list(map(lambda p: to_homogeneous(ref + to_cartesian(p), p[-1]), npw)))
        return cls(ref, to_cartesian(npw[0]), 360, norm_vec)

    @property
    def length(self):
        return self.radius * math.radians(self.theta)

    def curvature(self, u):
        return 1.0 / self.radius


class ConicArc(Crv):
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

        '''Setup'''
        super(ConicArc, self).__init__(nu, npw)


class LocalCubicInterpolatedCrv(Crv):
    def __init__(self, pts, tv):
        """
        Interpolate points using cubic bezier curve with specified tangent vector on each segment ending.
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


def point_inverse(c, p, dim=None, e1=1e-7):
    """
    Find the parameter 'u' s.t. c(u) = p
    :param c: Target curve.
    :type c: Crv
    :param p: Target point.
    :param dim: Dimension indicator.
    :type dim: int
    :param e1: Default error criteria.
    :type e1: float
    :return: The parameter.
    :rtype: float
    """

    if dim is not None and dim >= 3:
        raise ValueError("Inconsistent input.")

    '''Find initial u0'''
    val = np.unique(c.U)
    seg = len(val) - 1
    uc = []
    for i in range(seg):
        cs = val[i]
        ce = val[i + 1]
        cud = list(np.linspace(cs, ce, 10))
        uc += cud[:-1]
    uc.append(val[-1])

    min_idx = 0
    min_dist = sys.float_info.max

    if dim is None:
        for k, pu in enumerate(uc):
            cd = pnt_dist(c(pu), p)
            if cd < min_dist:
                min_dist = cd
                min_idx = k
        u0 = uc[min_idx]

        '''Newton Iteration'''
        tmp1 = c(u0) - p
        eps1 = norm(tmp1)
        tmp2 = c(u0, 1)
        eps2 = math.fabs(np.dot(tmp1, tmp2)) / (norm(tmp1) * norm(tmp2))
        while eps1 > e1 or eps2 > e1:
            u = u0 - np.dot(tmp2, tmp1) / (np.dot(c(u0, 2), tmp1) + norm(tmp2) ** 2)
            tmp1 = c(u) - p
            eps1 = norm(tmp1)
            tmp2 = c(u, 1)
            eps2 = math.fabs(np.dot(tmp1, tmp2)) / (norm(tmp1) * norm(tmp2))
            u0 = u
    else:
        for k, pu in enumerate(uc):
            cd = math.fabs(c(pu)[dim] - p)
            if cd < min_dist:
                min_dist = cd
                min_idx = k
        u0 = uc[min_idx]

        '''Newton Iteration'''
        tmp1 = c(u0)[dim] - p
        eps1 = math.fabs(tmp1)
        while eps1 > e1:
            u = u0 - tmp1 / c(u0, 1)[dim]
            tmp1 = c(u)[dim] - p
            eps1 = math.fabs(tmp1)
            u0 = u

    return u0


class NURBSCrvTester(unittest.TestCase):
    def test_construction(self):
        # knot, ctrl points, weights
        u_vec = [[0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1],
                 [0, 0, 0, 0.5, 1, 1, 1],
                 [0, 0, 0, 0, 1, 1, 1, 1],
                 [0, 0, 0, 0, 1, 1, 1, 1],
                 [0, 0, 1, 1]]
        pnt = [[[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1], [1, 0]],
               [[1, 0], [1, 1], [-1, 1], [-1, 0]],
               [[1, 0], [1, 2], [-1, 2], [-1, 0]],
               [[sqrt3 / 2, 1 / 2], [sqrt3, -3], [-sqrt3, -3], [-sqrt3 / 2, 1 / 2]],
               [[10, 10], [100, 100]]]
        w = [[1, 1 / sqrt2, 1, 1 / sqrt2, 1, 1 / sqrt2, 1, 1 / sqrt2, 1],
             [1, 0.5, 0.5, 1],
             [1, 1 / 3, 1 / 3, 1],
             [1, 1 / 6, 1 / 6, 1],
             [1, 1]]
        ans = ['circle1.igs', 'circle2.igs', 'circle3.igs', 'circle4.igs', 'line.igs']

        # Just used to avoid warning
        self.assertTrue(len(u_vec) == len(pnt) == len(w) == len(ans))

        iges_model = Model()
        for k in range(len(ans)):
            iges_model.clear()
            geom = Crv(u_vec[k], list(map(lambda _p, _w: to_homogeneous(np.append(_p, [0]), _w), pnt[k], w[k])))
            iges_model.add(geom.to_iges())
            iges_model.save(ans[k])
            print(repr(geom))

    def test_call(self):
        # knot, ctrl points, weights
        u_vec = [[0, 0, 0, 1, 2, 3, 3, 3],
                 [0, 0, 0, 1, 1, 1]]
        p = [[(0, 0, 0), (1, 1, 0), (3, 2, 0), (4, 1, 0), (5, -1, 0)],
             [(1, 0, 0), (1, 1, 0), (0, 1, 0)]]
        w = [[1, 4, 1, 1, 1],
             [1, 1, 2]]

        # u, derivative
        data = [[(0, 0), (1, 0), (3, 0)],
                [(0, 1), (0, 2), (1, 1)]]
        ans = [[(0, 0, 0), (7 / 5, 6 / 5, 0), (5, -1, 0)],
               [(0, 2, 0), (-4, 0, 0), (-1, 0, 0)]]

        # Just used to avoid warning
        self.assertTrue(len(u_vec) == len(p) == len(w) == len(data) == len(ans))

        for i in range(len(data)):
            crv = Crv(u_vec[i], list(map(lambda _p, _w: to_homogeneous(_p, _w), p[i], w[i])))
            cur_ans = list(map(lambda t: crv(t[0], t[1]), data[i]))
            np.testing.assert_array_equal(cur_ans, ans[i])

    def test_length(self):
        # knot, ctrl points, weights
        u_vec = [[0, 0, 3, 3],
                 [0, 0, 0, 1, 1, 1]]
        p = [[(0, 0, 0), (1, 1, 0)],
             [(1, 0, 0), (1, 1, 0), (0, 1, 0)]]
        w = [[1, 1],
             [1, 1, 2]]

        ans = [sqrt2, math.pi / 2]

        for i in range(len(ans)):
            crv = Crv(u_vec[i], list(map(lambda _p, _w: to_homogeneous(_p, _w), p[i], w[i])))
            cur_ans = crv.length
            self.assertTrue(math.isclose(cur_ans, ans[i]))

    def test_curvature(self):
        # knot, ctrl points, weights
        u_vec = [[0, 0, 3, 3],
                 [0, 0, 0, 1, 1, 1]]
        p = [[(0, 0, 0), (1, 1, 0)],
             [(1, 0, 0), (1, 1, 0), (0, 1, 0)]]
        w = [[1, 1],
             [1, 1, 2]]

        data = [[0, 1.5, 3],
                [0, 0.5, 1]]
        ans = [[0, 0, 0],
               [1, 1, 1]]

        for i in range(len(ans)):
            crv = Crv(u_vec[i], list(map(lambda _p, _w: to_homogeneous(_p, _w), p[i], w[i])))
            cur_ans = list(map(lambda u: crv.curvature(u), data[i]))
            for j in range(len(cur_ans)):
                self.assertTrue(math.isclose(cur_ans[j], ans[i][j]))

    def test_reverse(self):
        u_vec = [0, 0, 0, 1, 3, 6, 6, 8, 8, 8]
        s_vec = [0, 0, 0, 2, 2, 5, 7, 8, 8, 8]
        p = [(1, 2, 2), (2, 4, 8), (3, 9, 27), (4, 16, 64), (5, 25, 125), (6, 36, 216)]
        q = [(6, 36, 216), (5, 25, 125), (4, 16, 64), (3, 9, 27), (2, 4, 8), (1, 2, 2)]
        w = [1, 1, 1, 1, 1, 1]

        crv = Crv(u_vec, list(map(lambda _p, _w: to_homogeneous(_p, _w), p, w)))
        crv.reverse()
        self.assertTrue(np.array_equal(crv.U, s_vec))
        self.assertTrue(np.array_equal(crv.cpt, q))

    def test_pan(self):
        # knot, ctrl points, weights
        u_vec = [0, 0, 0, 1, 1, 1]
        p = [(1, 0, 0), (1, 1, 0), (0, 1, 0)]
        w = [1, 1, 2]

        pan_dir = (3, 4, 5)
        ans = [(4, 4, 5), (4, 5, 5), (3, 5, 5)]

        crv = Crv(u_vec, list(map(lambda _p, _w: to_homogeneous(_p, _w), p, w)))
        iges_model = Model()
        iges_model.add(crv.to_iges())
        crv.pan(pan_dir)
        iges_model.add(crv.to_iges())
        iges_model.save('test_pan.igs')

        self.assertTrue(np.array_equal(crv.cpt, ans))
        self.assertTrue(np.array_equal(crv.weight, w))

    def test_rotate(self):
        # knot, ctrl points, weights
        u_vec = [0, 0, 0, 1, 1, 1]
        p = [(1, 0, 0), (1, 1, 0), (0, 1, 0)]
        w = [1, 1, 2]

        # anchor point, rotation axis, rotation angle
        data = [[(0, 0, 0), (1, 1, 1), 30],
                [(0, 0, 0), (1, 1, 1), 45],
                [(0, 0, 0), (1, 1, 1), 90],
                [(0, 0, 0), (1, 1, 1), 180]]
        ans = [[(0.9106836, 1 / 3, -0.2440169), (2 / 3, 1.2440169, 0.0893164), (-0.2440169, 0.9106836, 1 / 3)],
               [(0.8047379, 0.5058793, -0.3106172), (0.4941207, 1.3106172, 0.1952621), (-0.3106172, 0.8047379, 0.5058793)],
               [(1 / 3, 0.9106836, -0.2440169), (0.0893164, 1.2440169, 2 / 3), (-0.2440169, 1 / 3, 0.9106836)],
               [(-1 / 3, 2 / 3, 2 / 3), (1 / 3, 1 / 3, 4 / 3), (2 / 3, -1 / 3, 2 / 3)]]

        # Just used to avoid warning
        self.assertTrue(len(data) == len(ans))

        for i, rt in enumerate(data):
            iges_model = Model()
            crv = Crv(u_vec, list(map(lambda _p, _w: to_homogeneous(_p, _w), p, w)))
            # print(crv)
            iges_model.add(crv.to_iges())
            crv.rotate(rt[0], rt[1], rt[2])
            # print(crv)
            iges_model.add(crv.to_iges())
            iges_model.save('test_rotate-{}.igs'.format(i))

            np.testing.assert_array_almost_equal(crv.cpt, ans[i])

    def test_insert_knot(self):
        u_vec = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        data = [(2.5, 1),  # 插入1个不在原节点矢量中的节点
                (2, 1),  # 插入1个在原节点矢量中的节点
                (2.5, 3),  # 插入1个不在原节点矢量中的节点3次
                (2, 2),  # 插入1个在原节点矢量中的节点2次
                (2.5, 4)]  # 插入1个不在原节点矢量中的节点4次
        ans = ['test_insert-0.igs',
               'test_insert-1.igs',
               'test_insert-2.igs',
               'test_insert-3.igs',
               None]

        for i in range(len(data)):
            crv = Crv(u_vec, pw)
            try:
                crv.insert_knot(data[i][0], data[i][1])
            except ValueError as e:
                print('Ok, illegal insertion detected with following msg:\n{}'.format(e))
            else:
                iges_model = Model()
                iges_model.add(crv.to_iges())
                iges_model.save(ans[i])

        # Should yield identical results
        crv1 = Crv(u_vec, pw)
        crv2 = Crv(u_vec, pw)
        crv1.insert_knot(2)
        crv1.insert_knot(2)
        crv2.insert_knot(2, 2)

        self.assertTrue(np.array_equal(crv1.U, crv2.U))
        self.assertTrue(np.array_equal(crv1.Pw, crv2.Pw))

    def test_refine(self):
        u_vec = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        data = [2.5, 2.5, 2.5]
        ans = [0, 0, 0, 0, 1, 2, 2.5, 2.5, 2.5, 3, 4, 5, 5, 5, 5]

        crv = Crv(u_vec, pw)
        iges_model = Model()
        iges_model.add(crv.to_iges())
        iges_model.save('test_refine_original.igs')
        crv.refine(data)
        iges_model.clear()
        iges_model.add(crv.to_iges())
        iges_model.save('test_refine_after.igs')
        self.assertTrue(np.array_equal(crv.U, ans))

    def test_elevate(self):
        # Knot , ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 1, 1, 1, 1]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        iges_model = Model()
        crv = Crv(u_vec, pw)
        iges_model.add(crv.to_iges())
        iges_model.save('test_elevate_before.igs')
        crv.elevate(2)
        iges_model.clear()
        iges_model.add(crv.to_iges())
        iges_model.save('test_elevate_after.igs')

        self.assertTrue(True)

    def test_reparameterization(self):
        # Knot, ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0]
        pw = [(0, 3.14, 0, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        param = [(2, 1, 3, 2),
                 (2, 0, 3, 1)]

        iges_model = Model()
        for k, dt in enumerate(param):
            iges_model.clear()
            crv = Crv(u_vec, pw)
            iges_model.add(crv.to_iges())
            crv.reparameterization(dt[0], dt[1], dt[2], dt[3])
            iges_model.add(crv.to_iges())
            iges_model.save('test_reparameterization-{}.igs'.format(k))

        self.assertTrue(True)

    def test_split(self):
        # Knot , ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 1, 1, 1, 1]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        iges_model = Model()
        crv = Crv(u_vec, pw)
        iges_model.add(crv.to_iges())
        iges_model.save('test_split_before.igs')

        iges_model.clear()
        crv_seg = Crv.split(crv, [0.2, 0.6])
        for seg in crv_seg:
            iges_model.add(seg.to_iges())
        iges_model.save('test_split_after.igs')

        self.assertTrue(True)

    def test_decompose(self):
        # Knot , ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 1, 1, 1, 1]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        iges_model = Model()
        crv = Crv(u_vec, pw)
        iges_model.add(crv.to_iges())
        iges_model.save('test_decompose_before.igs')

        iges_model.clear()
        crv_seg = Crv.decompose(crv)
        for seg in crv_seg:
            iges_model.add(seg.to_iges())
        iges_model.save('test_decompose_after.igs')

        self.assertTrue(True)

    def test_point_inverse(self):
        # Knot , ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 1, 1, 1, 1]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]
        crv = Crv(u_vec, pw)

        data = [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9, 1]
        for u in data:
            self.assertEqual(u, point_inverse(crv, crv(u)))

    def test_spline(self):
        # airfoil, degree, knot method
        data = [('M6', 3, 'chord'), ('M6', 3, 'centripetal'), ('M6', 5, 'chord'), ('M6', 5, 'centripetal'),
                ('NACA0012', 3, 'chord'), ('NACA0012', 3, 'centripetal'), ('NACA0012', 5, 'chord'), ('NACA0012', 5, 'centripetal'),
                ('RAE2822', 3, 'chord'), ('RAE2822', 3, 'centripetal'), ('RAE2822', 5, 'chord'), ('RAE2822', 5, 'centripetal')]

        iges_model = Model()
        for k, dt in enumerate(data):
            spline_foil = Spline(read_airfoil_pts(dt[0]))
            iges_model.clear()
            iges_model.add(spline_foil.to_iges())
            iges_model.save('test_spline-{}_{}_{}_{}.igs'.format(k, dt[0], dt[1], dt[2]))

        self.assertTrue(True)

    def test_global_interp(self):
        # airfoil, degree, knot method
        data = [('M6', 3, 'chord'), ('M6', 3, 'centripetal'), ('M6', 5, 'chord'), ('M6', 5, 'centripetal'),
                ('NACA0012', 3, 'chord'), ('NACA0012', 3, 'centripetal'), ('NACA0012', 5, 'chord'), ('NACA0012', 5, 'centripetal'),
                ('RAE2822', 3, 'chord'), ('RAE2822', 3, 'centripetal'), ('RAE2822', 5, 'chord'), ('RAE2822', 5, 'centripetal')]

        iges_model = Model()
        for k, dt in enumerate(data):
            crv = GlobalInterpolatedCrv(read_airfoil_pts(dt[0]), dt[1], dt[2])
            iges_model.clear()
            iges_model.add(crv.to_iges())
            iges_model.save('test_global_interp_crv-{}_{}_{}_{}.igs'.format(k, dt[0], dt[1], dt[2]))

        self.assertTrue(True)

    def test_local_interp(self):
        cr = 12
        spn = 21

        '''Front'''
        fsl = np.array([1.2, 2.3, 2.45, 0])
        fsl[-1] = spn - sum(fsl[:-1])
        alpha = np.radians([50, 45, 33, 28])
        fp = np.zeros((5, 3))
        for i in range(4):
            fp[i + 1] = fp[i]
            fp[i + 1][0] += fsl[i] * math.tan(alpha[i])
            fp[i + 1][2] += fsl[i]
        ftv = np.array([[0, 0, 1],
                        [math.sin(alpha[1]), 0, math.cos(alpha[1])],
                        [math.sin(alpha[1]), 0, math.cos(alpha[1])],
                        [math.sin(alpha[3]), 0, math.cos(alpha[3])],
                        [math.sin(alpha[3]), 0, math.cos(alpha[3])]])

        '''Tail'''
        tsl = np.array([1.6, 4.0, 1.8, 0])
        tsl[-1] = spn - sum(tsl[:-1])
        beta = np.radians([-15, -48, -20, 15])
        tp = np.zeros((5, 3))
        tp[0][0] = cr
        for i in range(4):
            tp[i + 1] = tp[i]
            tp[i + 1][0] += tsl[i] * math.tan(beta[i])
            tp[i + 1][2] += tsl[i]
        ttv = np.array([[0, 0, 1],
                        [math.sin(beta[1]), 0, math.cos(beta[1])],
                        [math.sin(beta[1]), 0, math.cos(beta[1])],
                        [math.sin(beta[3]), 0, math.cos(beta[3])],
                        [math.sin(beta[3]), 0, math.cos(beta[3])]])

        '''Display'''
        iges_model = Model()
        fc = LocalCubicInterpolatedCrv(fp, ftv)
        tc = LocalCubicInterpolatedCrv(tp, ttv)
        iges_model.add(fc.to_iges())
        iges_model.add(tc.to_iges())
        iges_model.save('lci.igs')

        self.assertTrue(True)

    def test_circle(self):
        # r, center, norm_vec
        data = [[0.3, (0, 0, 0), (0, 0, 1)],
                [0.3, (0, 0, 0), (0, 1, 0)],
                [0.3, (0, 0, 0), (1, 0, 0)],
                [0.3, (0, 0, 0), (1, 1, 1)],
                [1.5, (-1, -1, -1), (0, 0, 1)],
                [1.5, (-2.7, 0, 0), (0, 1, 0)],
                [1.5, (0, -3.14, 0), (1, 0, 0)],
                [1.5, (0, 0, -0.618), (1, 1, 1)]]

        iges_model = Model()
        for k, dt in enumerate(data):
            arc = Arc.closed_circle(dt[0], dt[1], dt[2])
            iges_model.clear()
            iges_model.add(arc.to_iges())
            iges_model.save('test_circle-{}_{}.igs'.format(k, dt[0]))

        self.assertTrue(True)

    def test_arc(self):
        self.assertTrue(True)

        # center, start, theta, norm_vec
        data = [[(0, 0, 0), (0.1, 0, 0), -20, (0, 0, 1)],
                [(0, 0, 0), (0.1, 0, 0), 0, (0, 0, 1)],
                [(0, 0, 0), (0.1, 0, 0), 3, (0, 0, 1)],
                [(0, 0, 0), (0.2, 0, 0), 45, (0, 0, 1)],
                [(0, 0, 0), (0.3, 0, 0), 65.2, (0, 0, 1)],
                [(0, 0, 0), (0.4, 0, 0), 90, (0, 0, 1)],
                [(0, 0, 0), (0.5, 0, 0), 120, (0, 0, 1)],
                [(0, 0, 0), (0.6, 0, 0), 135, (0, 0, 1)],
                [(0, 0, 0), (0.7, 0, 0), 150, (0, 0, 1)],
                [(0, 0, 0), (0.8, 0, 0), 180, (0, 0, 1)],
                [(0, 0, 0), (0.9, 0, 0), 195, (0, 0, 1)],
                [(0, 0, 0), (1.0, 0, 0), 225, (0, 0, 1)],
                [(0, 0, 0), (1.1, 0, 0), 240, (0, 0, 1)],
                [(0, 0, 0), (1.2, 0, 0), 270, (0, 0, 1)],
                [(0, 0, 0), (1.3, 0, 0), 315, (0, 0, 1)],
                [(0, 0, 0), (1.4, 0, 0), 324, (0, 0, 1)],
                [(0, 0, 0), (1.5, 0, 0), 360, (0, 0, 1)],
                [(0, 0, 0), (0.1, 0, 0), 1024, (0, 0, 1)],
                [(0.5, 1.2, 1.3), (0.1, 0, 0), 3, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.2, 0, 0), 45, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.3, 0, 0), 65.2, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.4, 0, 0), 90, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.5, 0, 0), 120, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.6, 0, 0), 135, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.7, 0, 0), 150, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.8, 0, 0), 180, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.9, 0, 0), 195, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.0, 0, 0), 225, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.1, 0, 0), 240, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.2, 0, 0), 270, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.3, 0, 0), 315, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.4, 0, 0), 324, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.5, 0, 0), 360, (1, 1, 1)]]

        iges_model = Model()
        for k, dt in enumerate(data):
            try:
                arc = Arc(dt[0], dt[1], dt[2], dt[3])
            except ValueError as e:
                print('Exception caught with msg: {}'.format(e))
            else:
                iges_model.add(arc.to_iges())

        iges_model.save('test_arc.igs')

    def test_arc_2pnt(self):
        self.assertTrue(True)

        # start, end, theta, norm_vec
        data = [[(0, 0, 500), (0, 0, 100), 180, (0, 1, 0)],
                [(0, 0, 50), (0, 0, 10), 180, (0, 1, 0)],
                [(0, 0, 5), (0, 0, 1), 180, (0, 1, 0)]]

        iges_model = Model()
        for k, dt in enumerate(data):
            arc = Arc.from_2pnt(dt[0], dt[1], dt[2], dt[3])
            iges_model.clear()
            iges_model.add(arc.to_iges())
            iges_model.save('test_arc_pnt-{}.igs'.format(k))

    def test_conic(self):
        self.assertTrue(True)

        # a,b of an ellipse
        data = [(10, 6), (20, 8), (50, 12), (10, 25)]

        iges_model = Model()
        for k, dt in enumerate(data):
            a, b = dt
            arc = ConicArc((0, -b, 0), (1, 0, 0), (a, 0, 0), (0, 1, 0), (a / sqrt2, -b / sqrt2, 0))  # 1/4 Ellipse Arc
            iges_model.add(arc.to_iges())
        iges_model.save('test_conic.igs')


if __name__ == '__main__':
    unittest.main()
