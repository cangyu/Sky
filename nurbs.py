import sys
import math
import unittest
import bisect
import numpy as np
from copy import deepcopy
from numpy.linalg import norm
from scipy.integrate import romberg
from scipy.interpolate import BSpline, make_interp_spline
from scipy.linalg import solve
from scipy.misc import comb
from iges import Model, Entity110, Entity126, Entity128
from transform import Quaternion, DCM
from misc import array_smart_copy, normalize, pnt_dist

sqrt2 = math.sqrt(2)
sqrt3 = math.sqrt(3)

"""
Implementation of the NURBS curve and surface utility.

Note:
All the NURBS notations are in the 'Clamped' format by default.

TODO:
(1) Optimize the calculation of derivatives inside Crv.__call__(u, d)
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


def find_span(n: int, p: int, u: float, u_vec):
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

    if u < u_vec[0] or u > u_vec[-1]:
        raise ValueError("Target parameter \'{}\' is not located within given knot vector.".format(u))

    left = p
    right = n + 1

    if math.isclose(u, u_vec[0]):
        return left
    elif math.isclose(u, u_vec[-1]):
        return right
    else:
        t = bisect.bisect_right(u_vec, u, left, right)
        return t if math.isclose(u, u_vec[t]) else t - 1


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

    if not np.array_equal(pans1, pans2):
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
        ans = [2, 3, 3, 5, 5, 6, 8, 8, 9, 9, 10]

        for i in range(len(data)):
            cur_ans = find_span(data[i][0], data[i][1], data[i][2], data[i][3])
            self.assertEqual(cur_ans, ans[i])

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
        :type crv: Crv
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
    return Crv(nU, nPw)


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

    def length(self):
        return pnt_dist(self.start, self.end)

    def curvature(self, u):
        return 0.0

    def to_iges(self):
        return Entity110(to_cartesian(self.Pw[0]), to_cartesian(self.Pw[-1]))


class Arc(Crv):
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

        '''Set-up'''
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


def point_inverse(c, p, dim=None):
    """
    Find the parameter 'u' s.t. c(u) = p
    :param c: Target curve.
    :type c: Crv
    :param p: Target point.
    :param dim: Dimension indicator.
    :type dim: int
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
        while eps1 > 1e-7 or eps2 > 1e-7:
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
        while eps1 > 1e-7:
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

        self.spl = []
        q = self.q
        for i in range(self.n + 1):
            self.spl.append(BSpline(self.V, self.Pw[i], q))

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
        for spl in self.spl:
            r.append(spl(v, l))

        rw = np.copy(r)
        spl = BSpline.construct_fast(self.U, rw, self.p)
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

    def __repr__(self):
        return '\nU Knot:\n{}\nV Knot:\n{}\nControl points:\n{}\n'.format(self.U, self.V, self.Pw)

    def __str__(self):
        ret = '\nClamped NURBS Surface\nDegree:({},{})'.format(self.p, self.q)
        ret += self.__repr__()
        return ret

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

    def mirror(self, axis):
        """
        Mirror the surface along specified axis.
        :param axis: Direction axis.
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
    def split(cls, surf, ubrk, vbrk):
        """
        将曲面分割成若干子部分
        :param surf: Surface to be split
        :type surf: ClampedNURBSSurf
        :param ubrk: breaking knot in u-direction
        :param vbrk: breaking knot in v-direction
        :return: Collection of split surf
        """

        cp = surf.p
        cq = surf.q

        '''Pre-check'''
        if len(ubrk) != 0 and (min(ubrk) <= 0 or max(ubrk) >= 1):
            raise AssertionError("Invalid input.")

        if len(vbrk) != 0 and (min(vbrk) <= 0 or max(vbrk) >= 1):
            raise AssertionError("Invalid input.")

        '''Copy back break knot info'''
        uspk = sorted(ubrk)
        uspk.append(1.0)
        vspk = sorted(vbrk)
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
                csrf = ClampedNURBSSurf(useg, vseg, cpt)
                csrf.standard_reparameterization()
                csrf_seg.append(csrf)
            ucpis += ucpn - 1
            ret.append(csrf_seg)

        return ret


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
        :type crv: Crv
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

        uknot = c1.U
        vknot = np.array([0, 0, 1, 1], float)

        '''Control points'''
        pw = np.zeros((len(c1.Pw), 2, 4))
        for i in range(len(c1.Pw)):
            pw[i][0] = np.copy(c1.Pw[i])
            pw[i][1] = np.copy(c2.Pw[i])

        super(RuledSurf, self).__init__(uknot, vknot, pw)


class RevolvedSurf(ClampedNURBSSurf):
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
        Pj = crv.cpt
        wj = crv.weight
        m = crv.n
        npw = np.zeros((len(u_knot) - 2 - 1, m + 1, 4), float)
        for j in range(m + 1):
            O = point_to_line(Pj[j], center, axis)
            X = Pj[j] - O
            r = norm(X)
            X = normalize(X)
            Y = np.cross(axis, X)
            npw[0][j] = crv.Pw[j]
            P0 = Pj[j]
            T0 = Y
            index = 0
            for i in range(1, narcs + 1):
                P2 = O + r * cosines[i] * X + r * sines[i] * Y
                npw[index + 2][j] = to_homogeneous(P2, wj[j])
                T2 = -sines[i] * X + cosines[i] * Y
                npw[index + 1][j] = to_homogeneous(line_intersection(P0, T0, P2, T2), wm * wj[j])
                index += 2
                if i < narcs:
                    P0 = P2
                    T0 = T2

        super(RevolvedSurf, self).__init__(u_knot, crv.U, npw)


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


class Skinned(ClampedNURBSSurf):
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
        Q = np.zeros((n + 1, m + 1, 3))
        Qw = np.zeros((n + 1, m + 1, 4))

        for i in range(n + 1):
            Q[i] = calc_ctrl_pts(v_knot, q, pnt[i], v_param[i])

        for i in range(n + 1):
            for j in range(m + 1):
                Qw[i][j] = to_homogeneous(Q[i][j])

        super(Skinned, self).__init__(u_knot, v_knot, Qw)


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


if __name__ == '__main__':
    unittest.main()
