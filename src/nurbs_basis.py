import unittest
import bisect
import math
import numpy as np
from misc import normalize

"""
Implementation of some basic utilities in NURBS.

Note:
All the NURBS notations are in the 'Clamped' format by default.
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


if __name__ == '__main__':
    unittest.main()
