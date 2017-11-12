import unittest
import math
import numpy as np
import sys
from copy import deepcopy
from numpy.linalg import norm
from scipy.integrate import romberg
from scipy.interpolate import BSpline, make_interp_spline
from scipy.linalg import solve
from scipy.misc import comb
from iges import Model, Entity110, Entity126, Entity128
from transform import Quaternion, DCM
from misc import angle_from_3pnt, array_smart_copy, normalize, pnt_dist, read_airfoil_pts, sqrt2, sqrt3

"""
Implementation of NURBS Basis, Crv, Surf.

Note:
All the NURBS notations are in the 'Clamped' format by default.
All the units follow SI by default.
"""


def to_homogeneous(pnt, w=1.0):
    """
    Convert cartesian coordinate into homogeneous format with specific weight.
    :param pnt: Original cartesian coordinate.
    :param w: Weight
    :type w: float
    :return: Homogeneous coordinate.
    """

    pw = np.zeros(len(pnt) + 1)
    pw[-1] = w

    if math.isclose(w, 0):
        for i in range(len(pnt)):
            pw[i] = pnt[i]
    else:
        for i in range(len(pnt)):
            pw[i] = pnt[i] * w

    return pw


def to_cartesian(pnt):
    """
    Convert homogeneous coordinates into cartesian format.
    The last component is considered as the weight.
    :param pnt: Original homogeneous coordinates.
    :return: Cartesian coordinates.
    """

    n = len(pnt) - 1
    p = np.zeros(n)

    if math.isclose(pnt[-1], 0.):
        for i in range(n):
            p[i] = pnt[i]
    else:
        for i in range(n):
            p[i] = pnt[i] / pnt[-1]

    return p


def line_intersection(p1, u1, p2, u2, with_ratio=False):
    """
    Calculate the intersection of 2 straight lines.
    s.t. p = p1 + alpha1 * u1 = p2 + alpha2 * u2
    :param p1: Point on the first line.
    :param u1: Direction vector of the first line.
    :param p2: Point on the second line.
    :param u2: Direction vector of the second line.
    :param with_ratio: To indicate if additional parameters are returned or not.
    :return: with_ratio ? (alpha1, alpha2, p) : p
    """

    cp1 = np.copy(p1)
    cp2 = np.copy(p2)
    cu1 = np.copy(u1)
    cu2 = np.copy(u2)
    dp = cp2 - cp1

    if len(cp1) != len(cu1) != len(cp2) != len(cu2):
        raise AssertionError("Inconsistent dimension!")
    if not cu1.any():
        raise AssertionError("Invalid U1 direction vector!")
    if not cu2.any():
        raise AssertionError("Invalid U2 direction vector!")
    if not np.cross(cu1, cu2).any():
        err_msg = "Two lines are parallel." if np.cross(dp, cu1).any() else "Two lines coincide with each other."
        raise AssertionError(err_msg)

    nu1 = normalize(cu1)
    nu2 = normalize(cu2)
    vu1 = cu2 - np.dot(nu1, cu2) * nu1
    vu2 = cu1 - np.dot(nu2, cu1) * nu2
    alpha1 = np.dot(dp, vu2) / np.dot(cu1, vu2)
    alpha2 = -np.dot(dp, vu1) / np.dot(cu2, vu1)
    pans1 = cp1 + alpha1 * cu1
    pans2 = cp2 + alpha2 * cu2

    if not math.isclose(np.linalg.norm(pans1 - pans2), 0, abs_tol=1e-12):
        raise AssertionError("No intersection.")

    return (alpha1, alpha2, pans1) if with_ratio else pans1


def point_to_line(target, center, axis):
    """
    Project a point onto a line.
    :param target: Point to be projected.
    :param center: A point on the line.
    :param axis: Direction vector of the line.
    :return: The projection point.
    """

    u = np.copy(normalize(axis))
    if not u.any():
        raise AssertionError('Invalid line direction vector.')

    t = np.copy(target)
    c = np.copy(center)
    return c + np.dot(t - c, u) * u


class Basis(object):
    def __init__(self, u, p):
        """
        BSpline basis functions.
        :param u: Knot vector.
        :param p: Degree of the basis functions.
        :type p: int
        """

        assert len(u) >= 2 * (p + 1)
        for k in range(p + 1):
            assert math.isclose(u[k], 0) and math.isclose(u[-k], 1)

        self.knot = np.copy(u)
        self.deg = p

    def __call__(self, *args, **kwargs):
        """
        Calculate value with given derivatives for all control pts or just a single pnt.
        :param args: u, d
        :param kwargs:
        :return: Value at different ctrl pts.
        """

        if len(args) == 0:
            raise ValueError('Too less arguments.')
        if len(args) == 1:
            u = args[0]
            if type(u) not in (int, float) or u < 0 or u > 1:
                raise TypeError('Wrong type for the first argument.')
        elif len(args) == 2:
            u = args[0]
            if type(u) not in (int, float) or u < 0 or u > 1:
                raise TypeError('Wrong input for \'u\'(1st).')
            d = args[1]
            if type(d) is not int or d < 0 or d > self.deg:
                raise ValueError('Wrong input for \'d\'(2nd).')
        else:
            raise ValueError('Too many arguments.')


def find_span(n, p, u, u_vec):
    """
    Determine the segment where parameter u is located.
    The binary-search algorithm is employed.
    :param n: The last index of the control-point sequence.
    :type n: int
    :param p: Degree of the basis function for B-Spline.
    :type p: int
    :param u: Target parameter.
    :type u: float
    :param u_vec: Knot vector.
    :return: The index s.t. u belongs to [u_vec[i], u_vec[i+1]).
    :rtype: int
    """

    if u == u_vec[n + 1]:  # Corner case: u=U[m]，将其节点区间的下标设为n
        return n

    low = p
    high = n + 1

    mid = int((high + low) / 2)
    while u < u_vec[mid] or u >= u_vec[mid + 1]:
        if u < u_vec[mid]:
            high = mid
        else:
            low = mid
        mid = int((high + low) / 2)

    return mid


def all_basis_val(u, p, u_vec):
    """
    Calculate all the value of basis function at given parameter(denoted as: bfv(i, p, u)).
    :param u: Target parameter.
    :type u: float
    :param p: Degree of the basis function for B-Spline.
    :type p: int
    :param u_vec: Knot vector.
    :return: n+1 elements, but only bfv(i-p,p,u)~bfv(i,p,u) are not zero.
    """

    m = len(u_vec) - 1
    n = m - p - 1
    i = find_span(n, p, u, u_vec)

    bfv = np.zeros(p + 1)
    left = np.zeros(p + 1)
    right = np.zeros(p + 1)

    bfv[0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - u_vec[i + 1 - j]
        right[j] = u_vec[i + j] - u
        saved = 0.0
        for r in range(j):
            temp = bfv[r] / (right[r + 1] + left[j - r])
            bfv[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        bfv[j] = saved

    ans = np.zeros(n + 1, float)
    for k in range(i - p, i + 1):
        ans[k] = bfv[k - (i - p)]

    return ans


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

        # TODO: Optimize the calculation of derivatives.

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


def calc_pnt_param(pts, method, with_dist=False):
    """
    计算每个插值点所对应的参数。
    :param pts: 插值点坐标序列
    :param method: 参数计算方法
    :type method: str
    :param with_dist: Indicate if the total distance is returned as well.
    :type with_dist: bool
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

    return (param, d) if with_dist else param


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
        for j in range(n + 1):
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


class Circle(Crv):
    def __init__(self, center, start, theta, norm_vec):
        """
        Spatial Circle.
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

        '''Basic Circle'''
        nu, npw = Circle.construct_planar(self.radius, self.theta)
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

        super(Circle, self).__init__(nu, npw)

    @classmethod
    def construct_planar(cls, radius, theta):
        """
        Circle on XY plane, with center at (0, 0, 0)
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
        nu, npw = Circle.construct_planar(radius, 360)
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
    def __init__(self, p0, t0, p2, t2, p):
        """
        Construct an open conic arc.
        :param p0: Starting point.
        :param t0: Tangent vector at start(Just used to indicate direction).
        :param p2: Ending point.
        :param t2: Tangent vector at end(Just used to indicate direction).
        :param p: A point on the arc.
        """

        p1, w1 = make_one_arc(p0, t0, p2, t2, p)

        if w1 <= -1.0:
            raise ValueError('Invalid arc.')
        ang = angle_from_3pnt(p0, p1, p2)
        if w1 >= 1.0:
            seg_num = 1
        else:
            if w1 > 0.0 and ang > 60:
                seg_num = 1
            elif w1 < 0.0 and ang > 90:
                seg_num = 4
            else:
                seg_num = 2

        n = 2 * seg_num
        j = n + 1
        u = np.empty(n + 4, float)
        for i in range(3):
            u[i] = 0.0
            u[j + i] = 1.0

        pw = np.empty((n + 1, 4), float)
        pw[0] = to_homogeneous(p0)
        pw[n] = to_homogeneous(p2)

        if seg_num == 1:
            pw[1] = to_homogeneous(p1, w1)
        else:
            q1, s, r1, wqr = split_arc(p0, p1, w1, p2)
            if seg_num == 2:
                pw[2] = to_homogeneous(s)
                pw[1] = to_homogeneous(q1, wqr)
                pw[3] = to_homogeneous(r1, wqr)
                u[3] = u[4] = 0.5
            else:
                pw[4] = to_homogeneous(s)
                w1 = wqr
                hq1, hs, hr1, wqr = split_arc(p0, q1, w1, s)
                pw[2] = to_homogeneous(hs)
                pw[1] = to_homogeneous(hq1, wqr)
                pw[3] = to_homogeneous(hr1, wqr)
                hq1, hs, hr1, wqr = split_arc(s, r1, w1, p2)
                pw[6] = to_homogeneous(hs)
                pw[5] = to_homogeneous(hq1, wqr)
                pw[7] = to_homogeneous(hr1, wqr)
                for i in range(2):
                    u[i + 3] = 0.25
                    u[i + 5] = 0.5
                    u[i + 7] = 0.75

        super(ConicArc, self).__init__(u, pw)


def make_one_arc(_p0, _t0, _p2, _t2, _p):
    """
    Construct a conic arc in one segment.
    Mainly, the job focus on calculating the middle point and corresponding weight.
    :param _p0: Starting point.
    :param _t0: Tangent vector at start(Just used to indicate direction).
    :param _p2: Ending point.
    :param _t2: Tangent vector at end(Just used to indicate direction).
    :param _p: A point on the arc.
    :return: The middle point and its weight.
    """

    p0 = np.copy(_p0)
    t0 = np.copy(_t0)
    p2 = np.copy(_p2)
    t2 = np.copy(_t2)
    p = np.copy(_p)
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

    return p1, w1


def split_arc(_p0, _p1, w1, _p2):
    """
    Split a 1-segment open conic arc into 2 parts.
    Ellipse with central angle equals 180 is treated especially, where w1==0.
    :param _p0: The 1st control point, its weight is 1 by default.
    :param _p1: The 2nd control point.
    :param w1: Weight of the 2nd control point.
    :type w1: float
    :param _p2: The 3rd control point, its weight is 1 by default.
    :return: new middle points on each part, the joint ending, and weight for the middle points.
    """

    p0 = np.copy(_p0)
    p1 = np.copy(_p1)
    p2 = np.copy(_p2)
    if math.isclose(w1, 0):
        s = p1
        tmp = (p0 - p2) / 2
        q1 = s + tmp
        r1 = s - tmp
        wqr = 1.0 / sqrt2
        return q1, s, r1, wqr
    else:
        t = 1 + w1
        q1 = (p0 + w1 * p1) / t
        r1 = (w1 * p1 + p2) / t
        wqr = math.sqrt(t / 2)
        s = (q1 + r1) / 2
        return q1, s, r1, wqr


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

    # TODO: Determine better converge criteria.

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


class Surf(object):
    def __init__(self, u, v, pts, rational=True):
        """
        Base class for NURBS Surface.
        :param u: u方向节点矢量, n+1个元素
        :param v: v方向节点矢量，m+1个元素
        :param pts: 齐次坐标序列，(n+1)x(m+1)个元素
        :param rational: Indicate if the control points are in homogeneous form.
        :type rational: bool
        """

        self.U = np.copy(u)
        self.V = np.copy(v)
        if rational:
            self.Pw = np.copy(pts)
        else:
            n, m, dim = pts.shape
            self.Pw = np.empty((n, m, 4))
            for i in range(n):
                for j in range(m):
                    self.Pw[i][j] = to_homogeneous(pts[i][j])

        self.spl = []
        for i in range(self.n + 1):
            self.spl.append(BSpline(self.V, self.Pw[i], self.q))

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
        :return: Degree in U direction.
        :rtype: int
        """

        return len(self.U) - self.n - 2

    @property
    def q(self):
        """
        :return: Degree in V direction.
        :rtype: int
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
        Control points.
        """

        ans = np.zeros((self.n + 1, self.m + 1, 3))
        for i in range(self.n + 1):
            for j in range(self.m + 1):
                ans[i][j] = to_cartesian(self.Pw[i][j])
        return ans

    def __repr__(self):
        return 'U Knot:\n{}\nV Knot:\n{}\nControl points:\n{}'.format(self.U, self.V, self.Pw)

    def __str__(self):
        ret = 'Clamped NURBS Surface\nDegree:({},{})\n'.format(self.p, self.q)
        ret += self.__repr__()
        return ret

    def __call__(self, u, v, k=0, l=0):
        """
        求在给定位置(u,v)处的导矢量
        :param u: U方向参数
        :param v: V方向参数
        :param k: U方向求导次数
        :param l: V方向求导次数
        :return: (u,v)处偏导矢量
        """

        # TODO: Implement BSpline basis function and replace scipy related module.

        r = []
        for spl in self.spl:
            r.append(spl(v, l))

        rw = np.copy(r)
        spl = BSpline.construct_fast(self.U, rw, self.p)
        pw = spl(u, k)
        return to_cartesian(pw)

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

        self.spl = []
        q = self.q
        for i in range(self.n + 1):
            self.spl.append(BSpline(self.V, self.Pw[i], q))

    def reverse(self, direction):
        """
        曲面反向
        """

        if direction not in ('U', 'V', 'UV'):
            raise ValueError('Invalid direction choice!')

        if direction in ('U', 'UV'):
            self.U = np.full(self.U.shape, self.U[0] + self.U[-1]) - self.U[::-1]
            self.Pw = self.Pw[::-1, :, :]

        if direction in ('V', 'UV'):
            self.V = np.full(self.V.shape, self.V[0] + self.V[-1]) - self.V[::-1]
            self.Pw = self.Pw[:, ::-1, :]

        self.reset(self.U, self.V, self.Pw)

    def swap(self):
        """
        交换UV方向节点矢量与控制点
        :return: None.
        """

        tmp = self.U[:]
        self.U = self.V[:]
        self.V = tmp
        self.Pw = np.transpose(self.Pw, (1, 0, 2))
        self.reset(self.U, self.V, self.Pw)

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

        self.reset(self.U, self.V, self.Pw)

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
        self.reset(self.U, self.V, self.Pw)

    def mirror(self, axis):
        """
        Mirror the surface along specified axis.
        :param axis: Direction axis.
        :type axis: str
        :return: None.
        """

        '''Defensive check'''
        if axis in ('X', 'x'):
            idx = 0
        elif axis in ('Y', 'y'):
            idx = 1
        elif axis in ('Z', 'z'):
            idx = 2
        else:
            raise ValueError("Invalid axis")

        '''Modify control points'''
        for i in range(self.n + 1):
            for j in range(self.m + 1):
                self.Pw[i][j][idx] *= -1

        '''Update'''
        self.reset(self.U, self.V, self.Pw)

    def to_iges(self, *args, **kwargs):
        """
        将曲面以IGES标准中第128号实体呈现
        :return: IGES representation of this surface.
        :rtype: Entity128
        """

        w = self.weight
        poly = 0 if (w != np.ones(w.shape)).any() else 1
        cpt = self.cpt

        form = kwargs['form'] if 'form' in kwargs else 0
        if len(args) == 4:
            closed_u, closed_v, periodic_u, periodic_v = args
        else:
            closed_u = closed_v = periodic_u = periodic_v = 0

        return Entity128(self.U, self.V, self.p, self.q, self.n, self.m, cpt, w, closed_u, closed_v, poly, periodic_u, periodic_v, self.U[0], self.U[-1], self.V[0], self.V[-1], form)

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
                cc = Crv(self.U, self.Pw[:, j])
                cc.insert_knot(uv, r)
                crv_list.append(cc)
                for i in range(self.n + 2):
                    npw[i][j] = np.copy(cc.Pw[i])

            self.reset(crv_list[0].U, self.V, npw)

        else:
            crv_list = []
            npw = np.zeros((self.n + 1, self.m + 2, 4))
            for i in range(self.n + 1):
                cc = Crv(self.V, self.Pw[i, :])
                cc.insert_knot(uv, r)
                crv_list.append(cc)
                for j in range(self.m + 2):
                    npw[i][j] = np.copy(cc.Pw[j])

            self.reset(self.U, crv_list[0].U, npw)

    def extract(self, direction, uv):
        """
        提取等参数线
        :param direction: 方向
        :param uv: 等参数值
        :return: 给定方向上的等参数线
        :rtype: Crv
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

            return Crv(self.V, nqw)

        else:
            npw = np.zeros((self.n + 1, 4))
            for i in range(self.n + 1):
                spl = BSpline(self.V, self.Pw[i, :, :], self.q)
                npw[i] = spl(uv)

            return Crv(self.U, npw)

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
                cc = Crv(self.U, self.Pw[:, j, :])
                cc.refine(extra_knot)
                crv_list.append(cc)
                for i in range(nh):
                    npw[i][j] = np.copy(cc.Pw[i])

            self.reset(crv_list[0].U, self.V, npw)

        else:
            mh = self.m + 1 + len(extra_knot)
            npw = np.zeros((self.n + 1, mh, 4))
            for i in range(self.n + 1):
                cc = Crv(self.V, self.Pw[i, :, :])
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
                cc = Crv(self.U, self.Pw[:, j, :])
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
                cc = Crv(self.V, self.Pw[i, :, :])
                cc.elevate(tv)
                crv_list.append(cc)

            mh = len(crv_list[0].Pw)
            npw = np.zeros((self.n + 1, mh, 4))
            for i in range(self.n + 1):
                for j in range(mh):
                    npw[i][j] = np.copy(crv_list[i].Pw[j])

            self.reset(self.U, crv_list[0].U, npw)

    def reparameterization(self, alpha, beta, gamma, delta, direction):
        """
        使用齐次有理函数将曲面重新参数化
             alpha * u + beta
        s = -------------------
             gamma * u + delta
        :param alpha: parameter
        :type alpha: float
        :param beta: parameter
        :type beta: float
        :param gamma: parameter
        :type gamma: float
        :param delta: parameter
        :type delta: float
        :param direction: parameter
        :type direction: str
        :return: None
        """

        if direction not in ('U', 'u', 'V', 'v'):
            raise AssertionError("Invalid direction parameter.")

        def g(x):
            return (alpha * x + beta) / (gamma * x + delta)

        def nbla(x):
            return gamma * x - alpha

        cpt = self.cpt
        npw = np.empty_like(self.Pw)
        wb = self.weight
        factor = 1.0

        if direction in ('U', 'u'):
            s = np.copy(list(map(g, self.U)))
            for k in range(self.p):
                factor *= nbla(s[k])
            for i in range(self.n + 1):
                factor /= nbla(s[i])
                factor *= nbla(s[i + self.p])
                for j in range(self.m + 1):
                    wb[i][j] *= factor
                    npw[i][j] = to_homogeneous(cpt[i][j], wb[i][j])
            self.reset(s, self.V, npw)
        else:
            t = np.copy(list(map(g, self.V)))
            for k in range(self.q):
                factor *= nbla(t[k])
            for j in range(self.m + 1):
                factor /= nbla(t[j])
                factor *= nbla(t[j + self.q])
                for i in range(self.n + 1):
                    wb[i][j] *= factor
                    npw[i][j] = to_homogeneous(cpt[i][j], wb[i][j])
            self.reset(self.U, t, npw)

    def standard_reparameterization(self):
        """
        将U,V两个方向的节点都统一到[0, 1]
        :return: None
        """

        '''U direction reparameterization'''
        a = self.U[0]
        b = self.U[-1]
        alpha = 1.0
        beta = -a
        gamma = 0.0
        delta = b - a
        self.reparameterization(alpha, beta, gamma, delta, 'U')

        '''V direction reparameterization'''
        a = self.V[0]
        b = self.V[-1]
        alpha = 1.0
        beta = -a
        gamma = 0.0
        delta = b - a
        self.reparameterization(alpha, beta, gamma, delta, 'V')

    @classmethod
    def split(cls, surf, u_brk, v_brk):
        """
        将曲面分割成若干子部分
        :param surf: Surface to be split
        :type surf: Surf
        :param u_brk: breaking knot in u-direction
        :param v_brk: breaking knot in v-direction
        :return: Collection of split surf
        """

        cp = surf.p
        cq = surf.q

        '''Pre-check'''
        if len(u_brk) != 0 and (min(u_brk) <= 0 or max(u_brk) >= 1):
            raise AssertionError("Invalid input.")

        if len(v_brk) != 0 and (min(v_brk) <= 0 or max(v_brk) >= 1):
            raise AssertionError("Invalid input.")

        '''Copy back break knot info'''
        uspk = sorted(u_brk)
        uspk.append(1.0)
        vspk = sorted(v_brk)
        vspk.append(1.0)

        '''Statistic current surf knot info'''
        uval, ucnt = np.unique(surf.U, return_counts=True)
        ukdt = dict(zip(uval, ucnt))

        vval, vcnt = np.unique(surf.V, return_counts=True)
        vkdt = dict(zip(vval, vcnt))

        '''Construct knot to be inserted'''
        uek = []
        for u in uspk:
            exist_cnt = ukdt.get(u) if u in ukdt else 0
            tc = cp - exist_cnt
            if tc > 0:
                for k in range(tc):
                    uek.append(u)

        vek = []
        for v in vspk:
            exist_cnt = vkdt.get(v) if v in vkdt else 0
            tc = cq - exist_cnt
            if tc > 0:
                for k in range(tc):
                    vek.append(v)

        '''Insert knots'''
        vsrf = deepcopy(surf)
        vsrf.refine('U', np.copy(uek))
        vsrf.refine('V', np.copy(vek))

        '''Build knot segment'''
        usdt = []
        uprev = 0.0
        ucki = cp + 1
        for u in uspk:
            cu = []
            for k in range(cp + 1):
                cu.append(uprev)
            while ucki < len(vsrf.U) and vsrf.U[ucki] <= u:
                cu.append(vsrf.U[ucki])
                ucki += 1
            if ucki < len(vsrf.U):
                cu.append(u)
            uprev = u
            usdt.append(cu)

        vsdt = []
        vprev = 0.0
        vcki = cq + 1
        for v in vspk:
            cv = []
            for k in range(cq + 1):
                cv.append(vprev)
            while vcki < len(vsrf.V) and vsrf.V[vcki] <= v:
                cv.append(vsrf.V[vcki])
                vcki += 1
            if vcki < len(vsrf.V):
                cv.append(v)
            vprev = v
            vsdt.append(cv)

        '''Extract control points'''
        ret = []
        ucpis = 0
        for useg in usdt:
            ucpn = len(useg) - cp - 1
            vcpis = 0
            csrf_seg = []
            for vseg in vsdt:
                vcpn = len(vseg) - cq - 1
                cpt = vsrf.Pw[ucpis:ucpis + ucpn, vcpis:vcpis + vcpn]
                vcpis += vcpn - 1
                csrf = Surf(useg, vseg, cpt)
                csrf.standard_reparameterization()
                csrf_seg.append(csrf)
            ucpis += ucpn - 1
            ret.append(csrf_seg)

        return ret


class GlobalInterpolatedSurf(Surf):
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

        u = np.zeros(n + 1)
        v = np.zeros(m + 1)
        u[-1] = 1.0
        v[-1] = 1.0

        '''Parameters of U direction'''
        dist = np.zeros((n + 1, m + 1))
        for j in range(m + 1):
            td = calc_pnt_param(pts[:, j], u_method)
            for i in range(n + 1):
                dist[i][j] = td[i]
        for i in range(n):
            u[i] = np.mean(dist[i])

        '''Parameters of V Direction'''
        for i in range(n + 1):
            td = calc_pnt_param(pts[i], v_method)
            for j in range(0, m + 1):
                dist[i][j] = td[j]
        for j in range(m):
            v[j] = np.mean(dist[:, j])

        '''Knot Vectors'''
        u_knot = calc_knot_vector(u, p)
        v_knot = calc_knot_vector(v, q)

        '''Control Points'''
        cr = np.zeros((n + 1, m + 1, dim))
        for j in range(m + 1):
            tp = calc_ctrl_pts(u_knot, p, pts[:, j], u)
            for i in range(n + 1):
                cr[i][j] = tp[i]

        cp = np.zeros((n + 1, m + 1, dim))
        for i in range(n + 1):
            cp[i] = calc_ctrl_pts(v_knot, q, cr[i], v)

        cpw = np.zeros((n + 1, m + 1, dim + 1))
        for i in range(n + 1):
            for j in range(m + 1):
                cpw[i][j] = to_homogeneous(cp[i][j])

        super(GlobalInterpolatedSurf, self).__init__(u_knot, v_knot, cpw)


class BilinearSurf(Surf):
    def __init__(self, p):
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

        :param p:4个角点, 2x2
        """

        u_vec = np.array([0, 0, 1, 1], float)
        v_vec = np.array([0, 0, 1, 1], float)

        ul, vl, dim = p.shape
        assert ul == 2 and vl == 2

        pw = np.ones((ul, vl, 4), float)
        for i in range(ul):
            for j in range(vl):
                for d in range(dim):
                    pw[i][j][d] = p[i][j][d]

        super(BilinearSurf, self).__init__(u_vec, v_vec, pw)


class ExtrudedSurf(Surf):
    def __init__(self, crv, direction):
        """
        拉伸曲面
        :param crv: Curve to be extruded.
        :type crv: Crv
        :param direction: Direction vector.
        """

        u_knot = np.copy(crv.U)
        v_knot = np.array([0, 0, 1, 1], float)
        n = len(crv.cpt)
        cpw = np.zeros((n, 2, 4))
        for i in range(n):
            cpw[i][0] = cpw[i][1] = np.copy(crv.Pw[i])
            w_dir = to_homogeneous(direction, cpw[i][0][3])
            for d in range(3):
                cpw[i][1][d] += w_dir[d]

        super(ExtrudedSurf, self).__init__(u_knot, v_knot, cpw)


class RuledSurf(Surf):
    def __init__(self, _c1, _c2):
        """
        生成V方向的直纹面,即两条曲线之间的线性插值
        :param _c1: 第1条曲线
        :type _c1: Crv
        :param _c2: 第2条曲线
        :type _c2: ClampedNURBSCrv
        """

        '''Not change original curve'''
        c1 = deepcopy(_c1)
        c2 = deepcopy(_c2)

        '''Check'''
        if not math.isclose(c1.U[0], c2.U[0]):
            raise ValueError('Incompatible starting knot!')
        if not math.isclose(c1.U[-1], c2.U[-1]):
            raise ValueError('Incompatible ending knot!')

        '''Knot vector'''
        p = max(c1.p, c2.p)
        c1.elevate(p - c1.p)
        c2.elevate(p - c2.p)

        if len(c1.U) != len(c2.U) or not math.isclose(norm(c1.U - c2.U), 0):
            all_knot = merge_knot(c1.U, c2.U)
            x1 = different_knot(all_knot, c1.U)
            x2 = different_knot(all_knot, c2.U)
            c1.refine(x1)
            c2.refine(x2)

        u_knot = c1.U
        v_knot = np.array([0, 0, 1, 1], float)

        '''Control points'''
        pw = np.zeros((len(c1.Pw), 2, 4))
        for i in range(len(c1.Pw)):
            pw[i][0] = np.copy(c1.Pw[i])
            pw[i][1] = np.copy(c2.Pw[i])

        super(RuledSurf, self).__init__(u_knot, v_knot, pw)


class RevolvedSurf(Surf):
    def __init__(self, center, axis, theta, crv):
        """
        曲线绕过指定点的轴线旋转指定角度得到的曲面
        :param center: 旋转中心
        :param axis: 旋转轴，正方向按右手法则给定
        :param theta: 旋转角度
        :type theta: float
        :param crv: 母线
        :type crv: Crv
        """

        while theta <= 0:
            theta += 360
        while theta > 360:
            theta -= 360

        '''Basic variables'''
        narcs = int(math.ceil(theta / 90))
        theta = math.radians(theta)
        delta_theta = theta / narcs
        wm = math.cos(delta_theta / 2)

        '''U Direction knots'''
        u_knot = np.zeros(2 * narcs + 4)
        u_knot[-1] = u_knot[-2] = u_knot[-3] = 1.0
        delta_knot = 1.0 / narcs
        for i in range(1, narcs):
            cur_index = 1 + 2 * i
            u_knot[cur_index] = u_knot[cur_index + 1] = i * delta_knot

        '''Pre-compute sine and cosine stuff'''
        sines = np.zeros(narcs + 1, float)
        cosines = np.ones(narcs + 1, float)
        angle = 0.0
        for i in range(1, narcs + 1):
            angle += delta_theta
            cosines[i] = math.cos(angle)
            sines[i] = math.sin(angle)

        '''Compute weighted control points on each line'''
        pj = crv.cpt
        wj = crv.weight
        m = crv.n
        npw = np.zeros((len(u_knot) - 2 - 1, m + 1, 4), float)
        for j in range(m + 1):
            oc = point_to_line(pj[j], center, axis)
            xd = pj[j] - oc
            r = norm(xd)
            xd = normalize(xd)
            yr = np.cross(axis, xd)
            npw[0][j] = crv.Pw[j]
            p0 = pj[j]
            t0 = yr
            index = 0
            for i in range(1, narcs + 1):
                p2 = oc + r * cosines[i] * xd + r * sines[i] * yr
                npw[index + 2][j] = to_homogeneous(p2, wj[j])
                t2 = -sines[i] * xd + cosines[i] * yr
                npw[index + 1][j] = to_homogeneous(line_intersection(p0, t0, p2, t2), wm * wj[j])
                index += 2
                if i < narcs:
                    p0 = p2
                    t0 = t2

        super(RevolvedSurf, self).__init__(u_knot, crv.U, npw)


class Coons(Surf):
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
        :type c0u: ClampedNURBSCrv
        :param c1u:沿U方向第2条曲线
        :type c1u: ClampedNURBSCrv
        :param c0v:沿V方向第1条曲线
        :type c0v: ClampedNURBSCrv
        :param c1v:沿V方向第2条曲线
        :type c1v: Crv
        """

        '''Check 4 corners'''
        assert math.isclose(norm(c0u(0) - c0v(0)), 0.0)
        assert math.isclose(norm(c0u(1) - c1v(0)), 0.0)
        assert math.isclose(norm(c1v(1) - c1u(1)), 0.0)
        assert math.isclose(norm(c0v(1) - c1u(0)), 0.0)

        '''Corner points'''
        s = np.zeros((2, 2, 3))
        s[0][0] = np.copy(c0u(0))
        s[0][1] = np.copy(c0v(1))
        s[1][0] = np.copy(c1v(0))
        s[1][1] = np.copy(c1u(1))

        '''Base surf'''
        r1 = RuledSurf(c0u, c1u)
        r2 = RuledSurf(c0v, c1v)
        r2.swap()
        t = BilinearSurf(s)

        '''Elevate to same order'''
        pu = max(r1.p, r2.p, t.p)
        pv = max(r1.q, r2.q, t.q)
        r1.elevate(pu - r1.p, pv - r1.q)
        r2.elevate(pu - r2.p, pv - r2.q)
        t.elevate(pu - t.p, pv - t.q)

        '''Unify knot vector'''
        xu = merge_knot(merge_knot(r1.U, r2.U), t.U)
        xv = merge_knot(merge_knot(r1.V, r2.V), t.V)

        xr1u = different_knot(xu, r1.U)
        xr2u = different_knot(xu, r2.U)
        xtu = different_knot(xu, t.U)

        xr1v = different_knot(xv, r1.V)
        xr2v = different_knot(xv, r2.V)
        xtv = different_knot(xv, t.V)

        r1.refine('U', xr1u)
        r1.refine('V', xr1v)

        r2.refine('U', xr2u)
        r2.refine('V', xr2v)

        t.refine('U', xtu)
        t.refine('V', xtv)

        '''Calculate new control points'''
        npw = r1.Pw + r2.Pw - t.Pw

        super(Coons, self).__init__(xu, xv, npw)


class Skinned(Surf):
    def __init__(self, crv, p, q, v_method='chord'):
        """
        蒙皮曲面，非有理
        :param crv: 非有理曲线集合
        :param p: 目标曲面u方向次数(曲线方向)
        :param q: 目标曲面v方向次数(展向)
        :param v_method: v方向插值方法
        """

        '''Promote all curve to p order'''
        crv_list = []
        for c in crv:
            cc = deepcopy(c)
            cc.elevate(p - cc.p)
            crv_list.append(cc)

        '''Merge all knot vectors in U direction'''
        u_knot = np.copy(crv_list[0].U)
        for i in range(1, len(crv_list)):
            u_knot = merge_knot(u_knot, crv_list[i].U)

        '''Unify all curve knot vector'''
        for c in crv_list:
            xu = different_knot(u_knot, c.U)
            c.refine(xu)

        '''Knot vector in V direction'''
        n = len(u_knot) - p - 2
        m = len(crv_list) - 1
        pnt = np.zeros((n + 1, m + 1, 3))
        for j in range(m + 1):
            for i in range(n + 1):
                pnt[i][j] = to_cartesian(crv_list[j].Pw[i])

        v_param = np.zeros((n + 1, m + 1))
        vk = []
        for i in range(n + 1):
            v_param[i] = calc_pnt_param(pnt[i], v_method)
            vk.append(calc_knot_vector(v_param[i], q))

        v_knot = np.mean(vk, axis=0)

        '''Calculate control points'''
        cq = np.zeros((n + 1, m + 1, 3))
        cqw = np.zeros((n + 1, m + 1, 4))

        for i in range(n + 1):
            cq[i] = calc_ctrl_pts(v_knot, q, pnt[i], v_param[i])

        for i in range(n + 1):
            for j in range(m + 1):
                cqw[i][j] = to_homogeneous(cq[i][j])

        super(Skinned, self).__init__(u_knot, v_knot, cqw)


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


class LocalBiCubicInterpSurf(Surf):
    def __init__(self, pts):
        """
        Locally bi-cubic interpolated bezier surface.
        :param pts: (n+1) x (m+1) points to be interpolated.
        """

        assert len(pts.shape) == 3
        assert pts.shape[-1] == 3

        '''Basic Variables'''
        dim = pts.shape[-1]
        n = pts.shape[0] - 1
        m = pts.shape[1] - 1

        '''Knot vector'''
        u_len = np.zeros(m + 1)
        v_len = np.zeros(n + 1)
        u = np.zeros(n + 1)
        v = np.zeros(m + 1)
        u_knot = np.zeros(2 * n + 6)
        v_knot = np.zeros(2 * m + 6)
        for i in range(1, 5):
            u_knot[-i] = v_knot[-i] = 1

        dist = np.zeros((m + 1, n + 1))
        for j in range(m + 1):
            dist[j], u_len[j] = calc_pnt_param(pts[:, j], 'chord', True)
        for i in range(n + 1):
            u_knot[2 * i + 2] = u_knot[2 * i + 3] = u[i] = np.mean(dist[:, i])
        dist = np.zeros((n + 1, m + 1))
        for i in range(n + 1):
            dist[i], v_len[i] = calc_pnt_param(pts[i], 'chord', True)
        for j in range(m + 1):
            v_knot[2 * j + 2] = v_knot[2 * j + 3] = v[j] = np.mean(dist[:, j])

        '''Unit tangent vector'''
        u_tan = np.empty((n + 1, m + 1, 3))
        v_tan = np.empty((n + 1, m + 1, 3))

        for j in range(m + 1):
            q = np.empty((n + 4, 3))
            for i in range(1, n + 1):
                q[i + 1] = pts[i][j] - pts[i - 1][j]
            q[1] = 2 * q[2] - q[3]
            q[0] = 2 * q[1] - q[2]
            q[n + 2] = 2 * q[n + 1] - q[n]
            q[n + 3] = 2 * q[n + 2] - q[n + 1]
            alpha = np.empty(n + 1)
            for i in range(n + 1):
                t1 = norm(np.cross(q[i], q[i + 1]))
                t2 = t1 + norm(np.cross(q[i + 2], q[i + 3]))
                alpha[i] = 1 if math.isclose(t2, 0) else t1 / t2
            for i in range(n + 1):
                u_tan[i][j] = normalize((1 - alpha[i]) * q[i + 1] + alpha[i] * q[i + 2])

        for i in range(n + 1):
            q = np.empty((m + 4, 3))
            for j in range(1, m + 1):
                q[j + 1] = pts[i][j] - pts[i][j - 1]
            q[1] = 2 * q[2] - q[3]
            q[0] = 2 * q[1] - q[2]
            q[m + 2] = 2 * q[m + 1] - q[m]
            q[m + 3] = 2 * q[m + 2] - q[m + 1]
            alpha = np.empty(m + 1)
            for j in range(m + 1):
                t1 = norm(np.cross(q[j], q[j + 1]))
                t2 = t1 + norm(np.cross(q[j + 2], q[j + 3]))
                alpha[i] = 1 if math.isclose(t2, 0) else t1 / t2
            for j in range(m + 1):
                v_tan[i][j] = normalize((1 - alpha[j]) * q[j + 1] + alpha[j] * q[j + 2])

        '''Control points'''
        cpt = np.zeros((2 * n + 2, 2 * m + 2, dim))
        cpt[0][0] = pts[0][0]
        cpt[-1][0] = pts[-1][0]
        cpt[0][-1] = pts[0][-1]
        cpt[-1][-1] = pts[-1][-1]

        for j in range(m + 1):
            for i in range(n):
                a = u_len[j] * (u[i + 1] - u[i]) / 3
                cpt[2 * i + 1][j] = pts[i][j] + a * u_tan[i][j]
                cpt[2 * i + 2][j] = pts[i + 1][j] - a * u_tan[i + 1][j]
        for i in range(n + 1):
            for j in range(m):
                a = v_len[i] * (v[j + 1] - v[j]) / 3
                cpt[i][2 * j + 1] = pts[i][j] + a * v_tan[i][j]
                cpt[i][2 * j + 2] = pts[i][j + 1] - a * v_tan[i][j + 1]

        # TODO: Calculate internal 4 ctrl_pts in each patch.

        super(LocalBiCubicInterpSurf, self).__init__(u_knot, v_knot, cpt, rational=False)


class BasicUtilityTester(unittest.TestCase):
    def test_find_span(self):
        # n, p, u, u_vec
        data = [[2, 2, 0, (0, 0, 0, 1, 1, 1)],
                [2, 2, 1, (0, 0, 0, 1, 1, 1)],
                [9, 3, 0.0, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.1, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.2, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.3, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.5, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.6, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.7, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.8, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 1.0, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)]]
        ans = [2, 2, 3, 5, 5, 6, 8, 8, 9, 9, 9]

        for k, dt in enumerate(data):
            cur_ans = find_span(dt[0], dt[1], dt[2], dt[3])
            self.assertEqual(cur_ans, ans[k])

    def test_all_basis_val(self):
        # u, p, u_vec
        data = [[2.5, 0, (0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5)],
                [2.5, 1, (0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5)],
                [2.5, 2, (0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5)]]
        ans = [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
               [0, 0, 0, 0.5, 0.5, 0, 0, 0, 0],
               [0, 0, 1 / 8, 3 / 4, 1 / 8, 0, 0, 0]]

        for i in range(len(data)):
            cur_ans = all_basis_val(data[i][0], data[i][1], data[i][2])
            for j in range(len(ans[i])):
                self.assertTrue(math.isclose(cur_ans[j], ans[i][j]))

    def test_line_intersection(self):
        # p1, u1, p2, u2
        data = [[(0, 5, 0), (1, 0, 0), (0, 5, 0), (0, 0, 1)],
                [(0, 0, 0), (1, 1, 0), (5, 0, 0), (1, -1, 0)],
                [(0, 1, 0), (0, 0, 1), (0, 2, 0), (1, 0, 0)]]
        ans = [(0, 5, 0),
               (2.5, 2.5, 0),
               None]

        for i in range(len(data)):
            try:
                cur_ans = line_intersection(data[i][0], data[i][1], data[i][2], data[i][3])
            except AssertionError as e:
                print("Exception caught! with msg: \'{}\'".format(e))
            else:
                for j in range(len(ans[i])):
                    self.assertTrue(math.isclose(cur_ans[j], ans[i][j]))

    def test_point_to_line(self):
        # t, p, u
        data = [[(0, 0, 0), (5, 0, 0), (1, -1, 0)],
                [(0, 2, 0), (0, 1, 0), (0, 0, 1)],
                [(3, 4, 5), (2, 2, 2), (0, 0, 0)]]

        ans = [(2.5, 2.5, 0),
               (0, 1, 0),
               None]

        for i in range(len(ans)):
            try:
                cur_ans = point_to_line(data[i][0], data[i][1], data[i][2])
            except AssertionError as e:
                print("Exception caught! with msg: \'{}\'".format(e))
            else:
                for j in range(len(ans[i])):
                    self.assertTrue(math.isclose(cur_ans[j], ans[i][j]))


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
            arc = Circle.closed_circle(dt[0], dt[1], dt[2])
            iges_model.clear()
            iges_model.add(arc.to_iges())
            iges_model.save('test_circle-{}_{}.igs'.format(k, dt[0]))

        self.assertTrue(True)

    def test_arc(self):
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
                arc = Circle(dt[0], dt[1], dt[2], dt[3])
            except ValueError as e:
                print('Exception caught with msg: {}'.format(e))
            else:
                iges_model.add(arc.to_iges())

        iges_model.save('test_arc.igs')
        self.assertTrue(True)

    def test_arc_2pnt(self):
        # start, end, theta, norm_vec
        data = [[(0, 0, 500), (0, 0, 100), 180, (0, 1, 0)],
                [(0, 0, 50), (0, 0, 10), 180, (0, 1, 0)],
                [(0, 0, 5), (0, 0, 1), 180, (0, 1, 0)]]

        iges_model = Model()
        for k, dt in enumerate(data):
            arc = Circle.from_2pnt(dt[0], dt[1], dt[2], dt[3])
            iges_model.clear()
            iges_model.add(arc.to_iges())
            iges_model.save('test_arc_pnt-{}.igs'.format(k))
        self.assertTrue(True)

    def test_conic(self):
        # a,b of an ellipse
        data = [(10, 6), (20, 8), (50, 12), (10, 25)]

        for k, dt in enumerate(data):
            a, b = dt
            arc1 = ConicArc((0, -b, 0), (1, 0, 0), (a, 0, 0), (0, 1, 0), (a / sqrt2, -b / sqrt2, 0))  # 1/4 Ellipse Circle
            arc2 = ConicArc((0, b, 0), (-1, 0, 0), (0, -b, 0), (1, 0, 0), (-a, 0, 0))  # 1/2 Ellipse Circle
            iges_model = Model()
            iges_model.add(arc1.to_iges())
            iges_model.add(arc2.to_iges())
            print(arc2)
            iges_model.save('test_conic-{}.igs'.format(k))
        self.assertTrue(True)


class NURBSSurfTester(unittest.TestCase):
    def test_call(self):
        pass

    def test_bilinear(self):
        l = 10.0
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        n = len(p)
        iges_model = Model()
        for k in range(n):
            s = BilinearSurf(p[k])
            iges_model.add(s.to_iges())
        iges_model.save('test_bilinear.igs')
        self.assertTrue(True)

    def test_pan(self):
        l = 10.0
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])
        d = [(0, 0, 10), (0, 0, 10)]

        self.assertTrue(len(p) == len(d))
        n = len(p)

        iges_model = Model()
        for k in range(n):
            iges_model.clear()
            s = BilinearSurf(p[k])
            iges_model.add(s.to_iges())
            s.pan(d[k])
            iges_model.add(s.to_iges())
            iges_model.save('test_pan-{}.igs'.format(k))

    def test_rotate(self):
        l = 10.0
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])
        ref = [[l, l, l], [l, l, l]]
        axis = [[0, 0, 5], [0, 0, 5]]
        ang = [45, 45]

        n = len(p)
        iges_model = Model()
        for k in range(n):
            iges_model.clear()
            s = BilinearSurf(p[k])
            iges_model.add(s.to_iges())
            s.rotate(ref[k], axis[k], ang[k])
            iges_model.add(s.to_iges())
            iges_model.save('test_rotate-{}.igs'.format(k))
        self.assertTrue(True)

    def test_reverse(self):
        l = 12.34
        p = np.array([[[[0, 0, 0], [0, l, l]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]],
                       [[l, 0, l], [l, l, 0]]]])
        axis = ['U', 'V', 'UV', 'U', 'V', 'UV']

        iges_model = Model()
        for k, bp in enumerate(p):
            iges_model.clear()
            s = BilinearSurf(bp)
            iges_model.add(s.to_iges())
            s.reverse(axis[k])
            iges_model.add(s.to_iges())
            iges_model.save('test_reverse-{}.igs'.format(k))
        self.assertTrue(True)

    def test_swap(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        for k, bpt in enumerate(p):
            s = BilinearSurf(bpt)
            ss = deepcopy(s)
            s.elevate(1, 2)
            ss.elevate(1, 2)
            s.swap()
            iges_model.clear()
            iges_model.add(s.to_iges())
            iges_model.add(ss.to_iges())
            iges_model.save('test_swap-{}.igs'.format(k))

            ni, nj = 50, 30
            u_dist, v_dist = np.meshgrid(np.linspace(0, 1, ni), np.linspace(0, 1, nj), indexing='ij')
            ps = np.zeros((ni, nj, 3))
            pss = np.zeros((ni, nj, 3))

            for i in range(ni):
                for j in range(nj):
                    ps[i][j] = s(v_dist[i][j], u_dist[i][j])
                    pss[i][j] = ss(u_dist[i][j], v_dist[i][j])

            self.assertTrue(np.allclose(ps, pss))

    def test_insert(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        print('Test insert utility.')
        for k, dt in enumerate(p):
            print('Test-{}'.format(k))
            iges_model.clear()
            s = BilinearSurf(dt)
            print('Before:{}'.format(s))
            iges_model.add(s.to_iges())
            s.insert_knot(0.5, 1, 'U')
            s.insert_knot(0.5, 1, 'V')
            print('After:{}'.format(s))
            iges_model.add(s.to_iges())
            iges_model.save('test_insert-{}.igs'.format(k))
        self.assertTrue(True)

    def test_refine(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        print('Test refine utility.')
        for k, dt in enumerate(p):
            print('Test-{}'.format(k))
            iges_model.clear()
            s = BilinearSurf(dt)
            ss = deepcopy(s)
            print('Before:{}'.format(s))
            iges_model.add(s.to_iges())
            s.refine('U', [0.2, 0.3, 0.6, 0.7])
            s.refine('V', [0.1, 0.4, 0.8, 0.9])
            print('After:{}'.format(s))
            iges_model.add(s.to_iges())
            iges_model.save('test_refine-{}.igs'.format(k))

            ni, nj = 50, 30
            u_dist = np.linspace(0, 1, ni)
            v_dist = np.linspace(0, 1, nj)
            ps = np.zeros((ni, nj, 3))
            pss = np.zeros((ni, nj, 3))

            for i in range(ni):
                for j in range(nj):
                    ps[i][j] = s(u_dist[i], v_dist[j])
                    pss[i][j] = ss(u_dist[i], v_dist[j])

            self.assertTrue(np.allclose(ps, pss))

    def test_elevate(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        print('Test elevate utility.')
        for k, dt in enumerate(p):
            print('Test-{}'.format(k))
            iges_model.clear()
            s = BilinearSurf(dt)
            ss = deepcopy(s)
            print('Before:{}'.format(s))
            iges_model.add(s.to_iges())
            s.elevate(1, 2)
            print('After:{}'.format(s))
            iges_model.add(s.to_iges())
            iges_model.save('test_elevate-{}.igs'.format(k))

            ni, nj = 50, 30
            u_dist = np.linspace(0, 1, ni)
            v_dist = np.linspace(0, 1, nj)
            ps = np.zeros((ni, nj, 3))
            pss = np.zeros((ni, nj, 3))
            for i in range(ni):
                for j in range(nj):
                    ps[i][j] = s(u_dist[i], v_dist[j])
                    pss[i][j] = ss(u_dist[i], v_dist[j])

            self.assertTrue(np.allclose(ps, pss))

    def test_split(self):
        print('Test split utility.')
        iges_model = Model()
        s = Surf(np.array([0, 0, 0, 0, 1., 1., 1., 1.]),
                 np.array([0, 0, 0, 0.5, 1, 1, 1]),
                 np.array([[[0, 0, 0, 1.], [0, 1, 1, 1], [0, 2, 3, 1], [0, 3, 2, 1]],
                           [[1, 0, 0, 1.], [1, 1, 2, 1], [1, 3, 5, 1], [1, 4, 2, 1]],
                           [[2, 0, 0, 1.], [2, 2, 2, 1], [2, 3, 6, 1], [2, 5, 7, 1]],
                           [[3, 0, 0, 1.], [3, 1, 1, 1], [3, 2, 2, 1], [3, 3, 3, 1]]]))
        iges_model.add(s.to_iges())
        print('Original:{}'.format(s))
        ss = Surf.split(s, (0.4, 0.5, 0.6), [0.2, 0.5, 0.7])
        k = 0
        for ue in ss:
            for ve in ue:
                iges_model.add(ve.to_iges())
                print('Segment {}:{}'.format(k, ve))
                k += 1
        iges_model.save('test_split.igs')
        self.assertTrue(True)

    def test_extract(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])
        ext = [('U', 0), ('U', 0.5), ('U', 1), ('V', 0), ('V', 0.5), ('V', 1),
               ('U', 0), ('U', 0.5), ('U', 1), ('V', 0), ('V', 0.5), ('V', 1)]
        self.assertTrue(len(p) == len(ext))

        iges_model = Model()
        print('Test extract utility.')
        for k, dt in enumerate(p):
            print('Case-{}'.format(k))
            s = BilinearSurf(dt)
            print('Surf:\n{}'.format(s))
            mid = s.extract(ext[k][0], ext[k][1])
            print('Curve at {}={}:\n{}'.format(ext[k][0], ext[k][1], mid))
            iges_model.clear()
            iges_model.add(s.to_iges())
            iges_model.add(mid.to_iges())
            iges_model.save('test_extract-{}_{}={}.igs'.format(k, ext[k][0], ext[k][1]))

    def test_extrude(self):
        crv = [Circle.from_2pnt((30, 45, 69), (44, 66, 88), 75, (1, -1, 1)),
               Circle.from_2pnt((0, 50, 0), (50, 0, 0), 315, (0, 0, 1))]
        direction = [(30, 30, 30), (0, 0, 90)]
        self.assertTrue(len(crv) == len(direction))

        iges_model = Model()
        for k in range(len(crv)):
            sf = ExtrudedSurf(crv[k], direction[k])
            iges_model.clear()
            iges_model.add(sf.to_iges())
            iges_model.save('test_extrude-{}.igs'.format(k))

    def test_revolved(self):
        generatrix = [Line((50, 0, 0), (50, 0, 20)),
                      Circle.from_2pnt((100, 0, 0), (100, 0, 20), 180, (1, 0, 0))]
        center = [(20, 0, 0), (50, 0, 0)]
        axis = [(0, 0, 1), (0, 0, 1)]
        theta = [60, 180]

        iges_model = Model()
        for k in range(len(generatrix)):
            iges_model.clear()
            iges_model.add(generatrix[k].to_iges())
            sf = RevolvedSurf(center[k], axis[k], theta[k], generatrix[k])
            iges_model.add(sf.to_iges())
            iges_model.save('test_revolved-{}.igs'.format(k))
        self.assertTrue(True)

    def test_ruled(self):
        # foil, z
        data = [('M6', 0), ('M6', 0.1), ('M6', 0.2), ('M6', 0.3),
                ('NACA0012', 0.4), ('NACA0012', 0.5), ('NACA0012', 0.6), ('NACA0012', 0.7),
                ('RAE2822', 0.8), ('RAE2822', 0.9), ('RAE2822', 1.0), ('RAE2822', 1.1)]
        crv = []
        for k in range(len(data)):
            foil = data[k][0]
            z = data[k][1]
            pts = read_airfoil_pts(foil)
            c = Spline(pts)
            c.pan((0, 0, z))
            crv.append(c)
        self.assertTrue(len(crv) == len(data))
        iges_model = Model()
        for i in range(1, len(crv)):
            sf = RuledSurf(crv[i - 1], crv[i])
            iges_model.add(sf.to_iges())
        iges_model.save('test_ruled.igs')

    def test_global_interp(self):
        # foil, N: num of profile, L: length per profile, p, q
        data = [('M6', 10, 0.7, 3, 3),
                ('M6', 10, 0.7, 3, 5),
                ('M6', 10, 0.7, 5, 3),
                ('M6', 10, 0.7, 5, 5),
                ('NACA0012', 10, 0.7, 3, 3),
                ('NACA0012', 10, 0.7, 3, 5),
                ('NACA0012', 10, 0.7, 5, 3),
                ('NACA0012', 10, 0.7, 5, 5),
                ('RAE2822', 10, 0.7, 3, 3),
                ('RAE2822', 10, 0.7, 3, 5),
                ('RAE2822', 10, 0.7, 5, 3),
                ('RAE2822', 10, 0.7, 5, 5)]

        iges_model = Model()
        for k in range(len(data)):
            foil = data[k][0]
            n = data[k][1]
            l = data[k][2]
            p = data[k][3]
            q = data[k][4]
            pts = read_airfoil_pts(foil)
            m, dim = pts.shape
            all_pts = np.zeros((n + 1, m, dim))
            for i in range(n + 1):
                all_pts[i] = np.copy(pts)
                for j in range(m):
                    all_pts[i][j][-1] = l * i
            wsf = GlobalInterpolatedSurf(all_pts, p, q)
            iges_model.clear()
            iges_model.add(wsf.to_iges())
            iges_model.save("test_global_interp-{}_{}_{}_{}.igs".format(k, foil, p, q))
        self.assertTrue(True)

    def test_local_interp(self):
        pass

    def test_coons(self):
        l = 20
        u0 = Spline(read_airfoil_pts('M6') * l)
        u0.pan((0, 0, 5))
        u1 = Circle.from_2pnt((25, 60, 5), (25, -60, 5), 180, (0, 0, 1))
        v0 = Line(u0.start, u1.start)
        v1 = Line(u0.end, u1.end)
        s = Coons(u0, u1, v0, v1)

        iges_model = Model()
        iges_model.add(s.to_iges())
        iges_model.add(u0.to_iges())
        iges_model.add(u1.to_iges())
        iges_model.add(v0.to_iges())
        iges_model.add(v1.to_iges())
        iges_model.save('test_coons.igs')
        self.assertTrue(True)

    def test_skinned(self):
        # foil, N: num of profile, L: length per profile, p, q
        data = [('M6', 10, 0.7, 3, 3),
                ('M6', 10, 0.7, 3, 5),
                ('M6', 10, 0.7, 5, 3),
                ('M6', 10, 0.7, 5, 5),
                ('NACA0012', 10, 0.7, 3, 3),
                ('NACA0012', 10, 0.7, 3, 5),
                ('NACA0012', 10, 0.7, 5, 3),
                ('NACA0012', 10, 0.7, 5, 5),
                ('RAE2822', 10, 0.7, 3, 3),
                ('RAE2822', 10, 0.7, 3, 5),
                ('RAE2822', 10, 0.7, 5, 3),
                ('RAE2822', 10, 0.7, 5, 5)]

        iges_model = Model()
        for k in range(len(data)):
            foil = data[k][0]
            n = data[k][1]
            l = data[k][2]
            p = data[k][3]
            q = data[k][4]
            pts = read_airfoil_pts(foil)
            m, dim = pts.shape
            all_pts = np.zeros((n + 1, m, dim))
            for i in range(n + 1):
                all_pts[i] = np.copy(pts)
                for j in range(m):
                    all_pts[i][j][-1] = l * i

            crv_list = [GlobalInterpolatedCrv(all_pts[i], p) for i in range(n + 1)]
            wsf = Skinned(crv_list, p, q)
            iges_model.clear()
            iges_model.add(wsf.to_iges())
            iges_model.save("test_Skinned-{}_{}_{}_{}.igs".format(k, foil, p, q))
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
