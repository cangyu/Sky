import os
import math
import numpy as np
from numpy.linalg import norm
from settings import AIRFOIL_DIR

sqrt2 = math.sqrt(2)
sqrt3 = math.sqrt(3)


def equal_check(*args):
    if len(args) < 2:
        return True

    prev = np.copy(args[0])
    for k in range(1, len(args)):
        crd = np.copy(args[k])
        if math.isclose(norm(prev - crd), 0.0, abs_tol=1e-8):
            prev = crd
        else:
            return False

    return True


def share(p, a, b):
    return (1 - p) * a + p * b


def angle_from_3pnt(p0, p1, p2):
    a = np.copy(p0) - np.copy(p1)
    b = np.copy(p2) - np.copy(p1)
    return math.degrees(math.acos(np.dot(a, b) / (norm(a) * norm(b))))


def vector_square(u):
    return sum(map(lambda x: x ** 2, u))


def array_smart_copy(src, dst):
    """
    Copy data without overflow.
    :param src: Source.
    :param dst: Destination.
    :return: None.
    """

    n = min(len(src), len(dst))
    for i in range(n):
        dst[i] = src[i]


def normalize(vec):
    """
    Normalize the input vector.
    :param vec: Original vector.
    :return: Normalized vector.
    """

    v = np.copy(vec)
    tmp = norm(v, 2)
    return v if math.isclose(tmp, 0) else v / tmp


def pnt_pan(origin, direction):
    return np.copy(origin) + np.copy(direction)


def pnt_dist(lhs, rhs):
    """
    Calculate the distance between two points with the great common degree.
    :param lhs: The 1st point.
    :param rhs: The 2nd point.
    :return: The distance between the 2 points.
    :rtype: float
    """

    ans = 0.
    dim = min(len(lhs), len(rhs))

    for i in range(dim):
        ans += math.pow(lhs[i] - rhs[i], 2)

    return math.sqrt(ans)


class XFoil(object):
    def __init__(self):
        pass
