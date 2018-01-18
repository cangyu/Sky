import math
import numpy as np
from scipy.optimize import newton

"""
Implementation of the Spacing control.

Note:
All the node distribute through [0, 1] by default.
"""


def uniform(*args):
    """
    Uniform distribution on [0,1]/[a,b].
    :param args: Num of nodes/(Starting, Ending, Num of nodes).
    :return: Uniform distribution on [0, 1] in numpy array form.
    """

    if len(args) == 1:
        n = args[0]
        return np.linspace(0.0, 1.0, n)
    else:
        a, b, n = args
        return np.linspace(a, b, n)


def chebshev_dist(start, end, n):
    """
    Simple chebshev distribution.
    :param start: Starting value.
    :type start: float
    :param end: Ending value.
    :type end: float
    :param n: Number of nodes.
    :type n: int
    :return: Chebshev distribution on [start, end] in numpy array form.
    """

    ang = np.linspace(np.pi, 0, n)
    pr = np.cos(ang)
    for i in range(0, n):
        pr[i] = start + (end - start) / 2 * (pr[i] + 1)

    return pr


def chebshev_dist_multi(seg, num):
    """
    Multi-segment Chebshev distribution.
    :param seg: Splitting values.
    :param num: Num of nodes within each segment.
    :return: Multi-segment Chebshev distribution in numpy array form.
    """

    assert len(seg) == len(num) + 1
    ret = chebshev_dist(seg[0], seg[1], num[0])
    for k in range(1, len(num)):
        csg = chebshev_dist(seg[k], seg[k + 1], num[k])
        ret = np.concatenate((ret, csg[1:]))

    return ret


def single_exponential(n, a):
    """
    Single exponential distribution of n nodes in [0, 1].
    :param n: Number of nodes.
    :type n: int
    :param a: Control parameter.
              >0 Nodes aggregate towards starting position.
              <0 Nodes aggregate towards ending position.
    :type a: float
    :return: Single exponential distribution of n nodes in [0, 1] in numpy array form.
    """

    assert not math.isclose(a, 0)

    t = (math.exp(a) - 1)
    return np.array([(math.exp(a * rho) - 1) / t for rho in uniform(n)])


def double_exponential(n, a1, a2, a3):
    """
    Double exponential distribution of n nodes in [0, 1].
    The curve pass through (a3, a1).
    :param n: Number of nodes.
    :type n: int
    :param a1: Horizontal coordinate of the control point.
    :type a1: float
    :param a2: Control parameter.
               >0 Nodes aggregate towards boundary points.
               <0 Nodes aggregate towards center.
    :type a2: float
    :param a3: Vertical coordinate of the control point.
    :type a3: float
    :return: Double exponential distribution of n nodes in [0, 1] in numpy array form.
    """

    p = (1 - a1) * a3 * (math.exp(a2) - 1) / ((1 - a3) * a1 * a2 * math.exp(a2))
    assert p > 0

    p1z = np.log(p)  # zero point of the 1st order derivative
    a4 = newton(func=lambda x: math.exp(x) - 1.0 - p * x, x0=5 * p1z, fprime=lambda x: math.exp(x) - p, maxiter=50, fprime2=lambda x: math.exp(x))

    ea21 = np.exp(a2) - 1
    ea41 = np.exp(a4) - 1

    def f(x):
        return a1 * (np.exp(a2 / a3 * x) - 1) / ea21 if x <= a3 else a1 + (1 - a1) * (np.exp(a4 / (1 - a3) * (x - a3)) - 1) / ea41

    return np.array([f(rho) for rho in uniform(n)])


def hyperbolic_tangent(n, b):
    """
    Hyperbolic tangent distribution of n nodes in [0, 1].
    :param n: Number of nodes.
    :type n: int
    :param b: Control parameter.
    :type b: float
    :return: Hyperbolic tangent distribution of n nodes in [0, 1] in numpy array form.
    """

    tb = math.tanh(b)

    def f(x):
        return 1 + math.tanh(b * (x - 1)) / tb

    return np.array([f(rho) for rho in uniform(n)])


def hyperbolic_sine(n, c):
    """
    Hyperbolic sine distribution of n nodes in [0, 1].
    :param n: Number of nodes.
    :type n: int
    :param c: Control parameter.
    :type c: float
    :return: Hyperbolic sine distribution of n nodes in [0, 1] in numpy array form.
    """

    sc = math.sinh(c)

    def f(x):
        return 1 + math.sinh(c * (x - 1)) / sc

    return np.array([f(rho) for rho in uniform(n)])


def linear_expand(seq, begin, end):
    """
    Rearrange the range of node sequence.
    :param seq: Node sequence to be rearranged.
    :param begin: New starting value.
    :type begin: float
    :param end: New ending value.
    :type end: float
    :return: Rearranged sequence.
    """

    r_min = min(seq)
    r_max = max(seq)
    dlt = end - begin
    r_dlt = r_max - r_min
    return np.full_like(seq, begin) if r_dlt == 0 else np.copy(list(map(lambda u: (u - r_min) / r_dlt * dlt + begin, seq)))
