import math
import numpy as np
from numpy.linalg import norm
import subprocess


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
