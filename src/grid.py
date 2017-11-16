import unittest
import sys
import platform
from abc import abstractmethod
from copy import deepcopy
from enum import Enum, unique
import math
import numpy as np
import scipy
from scipy.optimize import newton
from scipy.linalg import norm
from scipy import sparse
from scipy.sparse.linalg import dsolve
import matplotlib.pyplot as plt
from misc import share, equal_check, vector_square

"""
Implementation of the Spacing control.

Note:
All the node distribute through [0, 1] by default.
"""


def uniform(n):
    """
    Uniform distribution on [0,1].
    :param n: Number of nodes.
    :type n: int
    :return: Uniform distribution on [0, 1] in numpy array form.
    """

    return np.linspace(0.0, 1.0, n)


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


class SpacingTestCase(unittest.TestCase):
    def test_single_exponential(self):
        # n, A
        data = [(101, 5),
                (101, 3),
                (101, 1),
                (101, 0.05),
                (101, -0.05),
                (101, -1),
                (101, -3),
                (101, -5)]

        plt.figure()
        for i in range(len(data)):
            r = single_exponential(data[i][0], data[i][1])
            plt.plot(np.linspace(0.0, 1.0, data[i][0]), r, label='A={}'.format(data[i][1]))
        plt.xlabel('u')
        plt.ylabel('f(u)')
        plt.title('Single Exponential')
        plt.legend()
        plt.show()
        self.assertTrue(True)

    def test_double_exponential(self):
        # n, A1, A2, A3, A2取负值使得中间较密, 取正值使得中间稀疏，两边密集
        data = [(401, 0.5, -1.5, 0.5),
                (401, 0.5, 1.5, 0.5)]

        plt.figure()
        for i in range(len(data)):
            r = double_exponential(data[i][0], data[i][1], data[i][2], data[i][3])
            plt.plot(np.linspace(0.0, 1.0, data[i][0]), r, label='A1={}, A2={}, A3={}'.format(data[i][1], data[i][2], data[i][3]))
        plt.xlabel('u')
        plt.ylabel('f(u)')
        plt.title('Double Exponential')
        plt.legend()
        plt.show()
        self.assertTrue(True)

    def test_hyperbolic_tangent(self):
        # n, B
        data = [(201, 3), (201, 2), (201, 1)]

        plt.figure()
        for k, dt in enumerate(data):
            r = hyperbolic_tangent(dt[0], dt[1])
            plt.plot(np.linspace(0.0, 1.0, dt[0]), r, label='B={}'.format(dt[1]))
        plt.xlabel('u')
        plt.ylabel('f(u)')
        plt.title('Hyperbolic Tangent')
        plt.legend()
        plt.show()
        self.assertTrue(True)

    def test_hyperbolic_sine(self):
        # n, C
        data = [(201, 3), (201, 2), (201, 1)]

        plt.figure()
        for k, dt in enumerate(data):
            r = hyperbolic_sine(dt[0], dt[1])
            plt.plot(np.linspace(0.0, 1.0, dt[0]), r, label='C={}'.format(dt[1]))
        plt.xlabel('u')
        plt.ylabel('f(u)')
        plt.title('Hyperbolic Sine')
        plt.legend()
        plt.show()
        self.assertTrue(True)

    def test_newton_raphson(self):
        # pa1, pa2, pa3
        data = [(0.4, -1.2, 0.5)]
        ans = [0]

        self.assertTrue(len(data) == len(ans))

        for k, dt in enumerate(data):
            pa1, pa2, pa3 = dt
            p = (1 - pa1) * pa3 * (scipy.exp(pa2) - 1) / ((1 - pa3) * pa1 * pa2 * scipy.exp(pa2))
            p1z = math.log(p)

            def f(x):
                return scipy.exp(x) - 1.0 - p * x

            def pf(x):
                return scipy.exp(x) - p

            def ppf(x):
                return scipy.exp(x)

            fp1z = f(p1z)
            pa4 = newton(f, 5 * p1z, fprime=pf, maxiter=20, fprime2=ppf)
            print(p, p1z, fp1z, pa4, f(pa4))


"""
Implementation of the Trans-finite Interpolation(TFI) technique.

Note:
Only linear version are implemented.
"""


class TFI(object):
    def __init__(self):
        self.grid = None


class LinearTFI2D(TFI):
    def __init__(self, c1, c2, c3, c4):
        """
        2维无限插值(Linear)
        :param c1: 平行于x轴方向的第1条曲线，调用得到3维坐标点
        :param c2: 平行于y轴方向的第1条曲线，调用得到3维坐标点
        :param c3: 平行于x轴方向的第2条曲线，调用得到3维坐标点
        :param c4: 平行于y轴方向的第2条曲线，调用得到3维坐标点
        """

        super(LinearTFI2D, self).__init__()

        '''Defensive Check'''
        if not equal_check(c1(0), c2(0)):
            raise AssertionError("C1 start and C2 start not meet.")
        if not equal_check(c1(1), c4(0)):
            raise AssertionError("C1 end and C4 start not meet.")
        if not equal_check(c4(1), c3(1)):
            raise AssertionError("C4 end and C3 end not meet.")
        if not equal_check(c2(1), c3(0)):
            raise AssertionError("C2 end and C3 start not meet.")

        '''Copy back in case the parameters get changed outside'''
        self.c1 = deepcopy(c1)
        self.c2 = deepcopy(c2)
        self.c3 = deepcopy(c3)
        self.c4 = deepcopy(c4)

        '''Compute intersections in advance'''
        self.P12 = self.c1(0)  # c1与c2交点
        self.P14 = self.c1(1)  # c1与c4交点
        self.P23 = self.c3(0)  # c2与c3交点
        self.P34 = self.c3(1)  # c3与c4交点

    def U(self, u, v):
        return (1 - u) * self.c2(v) + u * self.c4(v)

    def V(self, u, v):
        return (1 - v) * self.c1(u) + v * self.c3(u)

    def UV(self, u, v):
        return (1 - u) * (1 - v) * self.P12 + u * v * self.P34 + (1 - u) * v * self.P23 + (1 - v) * u * self.P14

    def __call__(self, u, v):
        """
        曲面在(u,v)处的坐标
        :param u: U方向参数
        :type u: float
        :param v: V方向参数
        :type v: float
        :return: Coordinate at (u, v)
        """

        return self.U(u, v) + self.V(u, v) - self.UV(u, v)

    def calc_grid(self, pu, pv):
        """
        根据网格点的参数分布计算对应的坐标
        :param pu: 所有网格点的U方向参数值
        :param pv: 所有网格点的V方向参数值
        :return: None
        """

        '''Dimension on each axis'''
        d1 = len(pu)
        d2 = len(pv)

        '''Complement of each parameter distribution'''
        puc = np.ones(d1) - pu
        pvc = np.ones(d2) - pv

        '''Boundary coordinates'''
        pc1 = np.empty((d1, 3), float)
        pc3 = np.empty((d1, 3), float)
        pc2 = np.empty((d2, 3), float)
        pc4 = np.empty((d2, 3), float)

        '''Utilize pre-computed results'''
        pc1[0] = self.P12
        pc1[-1] = self.P14
        pc3[0] = self.P23
        pc3[-1] = self.P34
        pc2[0] = self.P12
        pc2[-1] = self.P23
        pc4[0] = self.P14
        pc4[-1] = self.P34

        '''Compute boundary coordinates in advance'''
        for i in range(1, d1 - 1):
            u = pu[i]
            pc1[i] = self.c1(u)
            pc3[i] = self.c3(u)

        for j in range(1, d2 - 1):
            v = pv[j]
            pc2[j] = self.c2(v)
            pc4[j] = self.c4(v)

        '''Transfinite interpolation'''
        self.grid = np.empty((d1, d2, 3), float)
        for i in range(d1):
            for j in range(d2):
                self.grid[i][j] = puc[i] * pc2[j] + pu[i] * pc4[j] + pvc[j] * pc1[i] + pv[j] * pc3[i]
                self.grid[i][j] -= puc[i] * pvc[j] * self.P12 + pu[i] * pv[j] * self.P34 + puc[i] * pv[j] * self.P23 + pvc[j] * pu[i] * self.P14


class LinearTFI3D(TFI):
    def __init__(self, s1, s2, s3, s4, s5, s6):
        """
        3维无限插值(Linear),右手系
        :param s1: 垂直于x轴的第1个曲面(Left)
        :param s2: 垂直于x轴的第2个曲面(Right)
        :param s3: 垂直于y轴的第1个曲面(Front)
        :param s4: 垂直于y轴的第2个曲面(Back)
        :param s5: 垂直于z轴的第1个曲面(Bottom)
        :param s6: 垂直于z轴的第2个曲面(Top)
        """

        super(LinearTFI3D, self).__init__()

        '''Defensive check'''
        assert equal_check(s1(0, 0), s3(0, 0), s5(0, 0))
        assert equal_check(s2(0, 0), s3(0, 1), s5(1, 0))
        assert equal_check(s2(1, 0), s4(0, 1), s5(1, 1))
        assert equal_check(s1(1, 0), s4(0, 0), s5(0, 1))
        assert equal_check(s1(0, 1), s3(1, 0), s6(0, 0))
        assert equal_check(s2(0, 1), s3(1, 1), s6(1, 0))
        assert equal_check(s2(1, 1), s4(1, 1), s6(1, 1))
        assert equal_check(s1(1, 1), s4(1, 0), s6(0, 1))

        self.s1 = deepcopy(s1)
        self.s2 = deepcopy(s2)
        self.s3 = deepcopy(s3)
        self.s4 = deepcopy(s4)
        self.s5 = deepcopy(s5)
        self.s6 = deepcopy(s6)

        self.c35 = lambda u: self.s5(u, 0)  # intersection of s5, s3
        self.c15 = lambda v: self.s5(0, v)  # intersection of s5, s1
        self.c45 = lambda u: self.s5(u, 1)  # intersection of s5, s4
        self.c25 = lambda v: self.s5(1, v)  # intersection of s5, s2
        self.c36 = lambda u: self.s6(u, 0)  # intersection of s6, s3
        self.c16 = lambda v: self.s6(0, v)  # intersection of s6, s1
        self.c46 = lambda u: self.s6(u, 1)  # intersection of s6, s4
        self.c26 = lambda v: self.s6(1, v)  # intersection of s6, s2
        self.c13 = lambda w: self.s3(w, 0)  # intersection of s1, s3
        self.c23 = lambda w: self.s2(0, w)  # intersection of s3, s2
        self.c24 = lambda w: self.s4(w, 1)  # intersection of s2, s4
        self.c14 = lambda w: self.s1(1, w)  # intersection of s4, s1

        self.p135 = self.c13(0)  # intersection of s1, s3, s5
        self.p235 = self.c23(0)  # intersection of s2, s3, s5
        self.p245 = self.c24(0)  # intersection of s2, s4, s5
        self.p145 = self.c14(0)  # intersection of s1, s4, s5
        self.p136 = self.c13(1)  # intersection of s1, s3, s6
        self.p236 = self.c23(1)  # intersection of s2, s3, s6
        self.p246 = self.c24(1)  # intersection of s2, s4, s6
        self.p146 = self.c14(1)  # intersection of s1, s4, s6

    def U(self, u, v, w):
        return (1 - u) * self.s1(v, w) + u * self.s2(v, w)

    def V(self, u, v, w):
        return (1 - v) * self.s3(w, u) + v * self.s4(w, u)

    def W(self, u, v, w):
        return (1 - w) * self.s5(u, v) + w * self.s6(u, v)

    def UV(self, u, v, w):
        return (1 - u) * (1 - v) * self.c13(w) + (1 - u) * v * self.c14(w) + u * (1 - v) * self.c23(w) + u * v * self.c24(w)

    def VW(self, u, v, w):
        return (1 - v) * (1 - w) * self.c35(u) + (1 - v) * w * self.c36(u) + v * (1 - w) * self.c45(u) + v * w * self.c46(u)

    def WU(self, u, v, w):
        return (1 - w) * (1 - u) * self.c15(v) + (1 - w) * u * self.c25(v) + w * (1 - u) * self.c16(v) + w * u * self.c26(v)

    def UVW(self, u, v, w):
        return (1 - u) * (1 - v) * (1 - w) * self.p135 + (1 - u) * (1 - v) * w * self.p136 + (1 - u) * v * (1 - w) * self.p145 + (1 - u) * v * w * self.p146 + u * (1 - v) * (1 - w) * self.p235 + u * (1 - v) * w * self.p236 + u * v * (1 - w) * self.p245 + u * v * w * self.p246

    @classmethod
    def from_edges(cls, *args):
        """
        Construct 3D TFI from its boundary edges.
        :param args: The set of boundary edges
        :return: 3D TFI
        :rtype: LinearTFI3D
        """

        assert len(args) == 12

        s1 = LinearTFI2D(args[3], args[8], args[7], args[11])
        s2 = LinearTFI2D(args[1], args[9], args[5], args[10])
        s3 = LinearTFI2D(args[8], args[0], args[9], args[4])
        s4 = LinearTFI2D(args[11], args[2], args[10], args[6])
        s5 = LinearTFI2D(args[0], args[3], args[2], args[1])
        s6 = LinearTFI2D(args[4], args[7], args[6], args[5])

        return cls(s1, s2, s3, s4, s5, s6)

    def __call__(self, u, v, w):
        """
        在(u,v,w)处的坐标
        """

        return self.U(u, v, w) + self.V(u, v, w) + self.W(u, v, w) - (self.UV(u, v, w) + self.VW(u, v, w) + self.WU(u, v, w)) + self.UVW(u, v, w)

    def calc_grid(self, pu, pv, pw):
        """
        根据网格点的参数分布计算对应的坐标
        :param pu: 所有网格点的U方向参数值
        :param pv: 所有网格点的V方向参数值
        :param pw: 所有网格点的W方向参数值
        :return: None
        """

        '''Dimension of each axis'''
        d1 = len(pu)
        d2 = len(pv)
        d3 = len(pw)

        '''Complement of each parameter distribution'''
        puc = np.ones(d1) - pu
        pvc = np.ones(d2) - pv
        pwc = np.ones(d3) - pw

        '''Boundary edge coordinates'''
        pc35 = np.empty((d1, 3), float)
        pc36 = np.empty((d1, 3), float)
        pc46 = np.empty((d1, 3), float)
        pc45 = np.empty((d1, 3), float)
        pc16 = np.empty((d2, 3), float)
        pc26 = np.empty((d2, 3), float)
        pc25 = np.empty((d2, 3), float)
        pc15 = np.empty((d2, 3), float)
        pc13 = np.empty((d3, 3), float)
        pc23 = np.empty((d3, 3), float)
        pc24 = np.empty((d3, 3), float)
        pc14 = np.empty((d3, 3), float)

        '''Boundary surface coordinates'''
        ps1 = np.empty((d2, d3, 3), float)
        ps2 = np.empty((d2, d3, 3), float)
        ps3 = np.empty((d3, d1, 3), float)
        ps4 = np.empty((d3, d1, 3), float)
        ps5 = np.empty((d1, d2, 3), float)
        ps6 = np.empty((d1, d2, 3), float)

        '''Utilize pre-computed corner coordinates'''
        pc35[0] = self.p135
        pc35[-1] = self.p235
        pc15[0] = self.p135
        pc15[-1] = self.p145
        pc45[0] = self.p145
        pc45[-1] = self.p245
        pc25[0] = self.p235
        pc25[-1] = self.p245
        pc36[0] = self.p136
        pc36[-1] = self.p236
        pc16[0] = self.p136
        pc16[-1] = self.p146
        pc46[0] = self.p146
        pc46[-1] = self.p246
        pc26[0] = self.p236
        pc26[-1] = self.p246
        pc13[0] = self.p135
        pc13[-1] = self.p136
        pc23[0] = self.p235
        pc23[-1] = self.p236
        pc24[0] = self.p245
        pc24[-1] = self.p246
        pc14[0] = self.p145
        pc14[-1] = self.p146

        '''Compute coordinates on each edge in advance'''
        for i in range(1, int(d1 - 1)):
            u = pu[i]
            pc35[i] = self.c35(u)
            pc36[i] = self.c36(u)
            pc46[i] = self.c46(u)
            pc45[i] = self.c45(u)

        for j in range(1, int(d2 - 1)):
            v = pv[j]
            pc16[j] = self.c16(v)
            pc26[j] = self.c26(v)
            pc25[j] = self.c25(v)
            pc15[j] = self.c15(v)

        for k in range(1, int(d3 - 1)):
            w = pw[k]
            pc13[k] = self.c13(w)
            pc23[k] = self.c23(w)
            pc24[k] = self.c24(w)
            pc14[k] = self.c14(w)

        '''Copy edge coordinates to surface coordinate array'''
        ps1[0] = pc13
        ps1[-1] = pc14
        ps2[0] = pc23
        ps2[-1] = pc24
        for j in range(1, int(d2 - 1)):
            ps1[j][0] = pc15[j]
            ps1[j][-1] = pc16[j]
            ps2[j][0] = pc25[j]
            ps2[j][-1] = pc26[j]

        ps3[0] = pc35
        ps3[-1] = pc36
        ps4[0] = pc45
        ps4[-1] = pc46
        for k in range(1, d3 - 1):
            ps3[k][0] = pc13[k]
            ps3[k][-1] = pc23[k]
            ps4[k][0] = pc14[k]
            ps4[k][-1] = pc24[k]

        ps5[0] = pc15
        ps5[-1] = pc25
        ps6[0] = pc16
        ps6[-1] = pc26
        for i in range(1, d1 - 1):
            ps5[i][0] = pc35[i]
            ps5[i][-1] = pc45[i]
            ps6[i][0] = pc36[i]
            ps6[i][-1] = pc46[i]

        '''Compute surface inner coordinates'''
        for j in range(1, int(d2 - 1)):
            v = pv[j]
            for k in range(1, int(d3 - 1)):
                w = pw[k]
                ps1[j][k] = self.s1(v, w)
                ps2[j][k] = self.s2(v, w)

        for k in range(1, int(d3 - 1)):
            w = pw[k]
            for i in range(1, int(d1 - 1)):
                u = pu[i]
                ps3[k][i] = self.s3(w, u)
                ps4[k][i] = self.s4(w, u)

        for i in range(1, int(d1 - 1)):
            u = pu[i]
            for j in range(1, int(d2 - 1)):
                v = pv[j]
                ps5[i][j] = self.s5(u, v)
                ps6[i][j] = self.s6(u, v)

        '''Compute the grid'''
        self.grid = np.empty((d1, d2, d3, 3), float)
        for i in range(d1):
            for j in range(d2):
                for k in range(d3):
                    self.grid[i][j][k] = puc[i] * ps1[j][k] + pu[i] * ps2[j][k]
                    self.grid[i][j][k] += pvc[j] * ps3[k][i] + pv[j] * ps4[k][i]
                    self.grid[i][j][k] += pwc[k] * ps5[i][j] + pw[k] * ps6[i][j]
                    self.grid[i][j][k] -= puc[i] * pvc[j] * pc13[k] + puc[i] * pv[j] * pc14[k] + pu[i] * pvc[j] * pc23[k] + pu[i] * pv[j] * pc24[k]
                    self.grid[i][j][k] -= pvc[j] * pwc[k] * pc35[i] + pvc[j] * pw[k] * pc36[i] + pv[j] * pwc[k] * pc45[i] + pv[j] * pw[k] * pc46[i]
                    self.grid[i][j][k] -= pwc[k] * puc[i] * pc15[j] + pwc[k] * pu[i] * pc25[j] + pw[k] * puc[i] * pc16[j] + pw[k] * pu[i] * pc26[j]
                    self.grid[i][j][k] += puc[i] * pvc[j] * pwc[k] * self.p135 + puc[i] * pvc[j] * pw[k] * self.p136 + puc[i] * pv[j] * pwc[k] * self.p145 + puc[i] * pv[j] * pw[k] * self.p146 + pu[i] * pvc[j] * pwc[k] * self.p235 + pu[i] * pvc[j] * pw[k] * self.p236 + pu[i] * pv[j] * pwc[k] * self.p245 + pu[i] * pv[
                        j] * pw[k] * self.p246


class TFITestCase(unittest.TestCase):
    def test_2d_rect(self):
        # L, W, U, V
        data = [(5, 4, 11, 9),
                (8, 8, 31, 21)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([dt[0] * u, 0, 0]),  # C1
                              lambda v: np.array([0, dt[1] * v, 0]),  # C3
                              lambda u: np.array([dt[0] * u, dt[1], 0]),  # C2
                              lambda v: np.array([dt[0], dt[1] * v, 0]))  # C4
            tfi.calc_grid(np.linspace(0, 1, dt[2]),
                          np.linspace(0, 1, dt[3]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_rect-{}.xyz'.format(k))

    def test_2d_circle(self):
        # R1, R2, U, V
        data = [(1, 2, 6, 11),
                (0, 5, 16, 33)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([(1 - u) * dt[0] + u * dt[1], 0, 0]),
                              lambda v: np.array([dt[0] * math.cos(0.5 * math.pi * v), dt[0] * math.sin(0.5 * math.pi * v), 0]),
                              lambda u: np.array([0, (1 - u) * dt[0] + u * dt[1], 0]),
                              lambda v: np.array([dt[1] * math.cos(0.5 * math.pi * v), dt[1] * math.sin(0.5 * math.pi * v), 0]))

            tfi.calc_grid(np.linspace(0, 1, dt[2]),
                          np.linspace(0, 1, dt[3]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_circle-{}.xyz'.format(k))

    def test_2d_eccentric(self):
        # Delta, R1, R2, U, V
        data = [(-10, 4, 25, 16, 41),
                (-10, 0, 25, 16, 41),
                (-1, 2, 5, 16, 33)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([(dt[0] + dt[1]) * (1 - u) + dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[1] * math.cos(math.pi * v) + dt[0], dt[1] * math.sin(math.pi * v), 0]),
                              lambda u: np.array([(dt[0] - dt[1]) * (1 - u) - dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[2] * math.cos(math.pi * v), dt[2] * math.sin(math.pi * v), 0]))

            tfi.calc_grid(np.linspace(0, 1, dt[3]),
                          np.linspace(0, 1, dt[4]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_eccentric-{}.xyz'.format(k))

    def test_2d_crv_rect(self):
        # L, H1, H2, H3
        data = [(100, 40, 60, 10, 50, 25)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([u * dt[0], 4 * dt[3] * u * (1 - u), 0]),
                              lambda v: np.array([0, v * dt[1], 0]),
                              lambda u: np.array([u * dt[0], (dt[1] * (1 - u * u) + dt[2] * u * u), 0]),
                              lambda v: np.array([dt[0], v * dt[2], 0]))

            tfi.calc_grid(np.linspace(0, 1, dt[4]),
                          np.linspace(0, 1, dt[5]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_crv_rect-{}.xyz'.format(k))

    def test_3d_cuboid(self):
        # L, W, H, U, V, W
        data = [(5, 4, 3, 11, 9, 5),
                (10, 10, 10, 21, 31, 41)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI3D(lambda v, w: np.array([0, v * dt[1], w * dt[2]]),
                              lambda v, w: np.array([dt[0], v * dt[1], w * dt[2]]),
                              lambda w, u: np.array([u * dt[0], 0, w * dt[2]]),
                              lambda w, u: np.array([u * dt[0], dt[1], w * dt[2]]),
                              lambda u, v: np.array([u * dt[0], v * dt[1], 0]),
                              lambda u, v: np.array([u * dt[0], v * dt[1], dt[2]]))
            tfi.calc_grid(np.linspace(0, 1.0, dt[3]),
                          np.linspace(0, 1.0, dt[4]),
                          np.linspace(0, 1.0, dt[5]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_3d_cuboid-{}.xyz'.format(k))

    def test_3d_sect(self):
        # R_MIN, R_MAX, THETA_MIN, THETA_MAX, H_MIN, H_MAX, U, V, W
        data = [(5, 20, 60, 120, -50, 50, 31, 16, 61)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI3D(lambda v, w: np.array([dt[0] * math.cos(math.radians(share(v, dt[2], dt[3]))), dt[0] * math.sin(math.radians(share(v, dt[2], dt[3]))), share(w, dt[4], dt[5])]),
                              lambda v, w: np.array([dt[1] * math.cos(math.radians(share(v, dt[2], dt[3]))), dt[1] * math.sin(math.radians(share(v, dt[2], dt[3]))), share(w, dt[4], dt[5])]),
                              lambda w, u: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(dt[2])), share(u, dt[0], dt[1]) * math.sin(math.radians(dt[2])), share(w, dt[4], dt[5])]),
                              lambda w, u: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(dt[3])), share(u, dt[0], dt[1]) * math.sin(math.radians(dt[3])), share(w, dt[4], dt[5])]),
                              lambda u, v: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(share(v, dt[2], dt[3]))), share(u, dt[0], dt[1]) * math.sin(math.radians(share(v, dt[2], dt[3]))), dt[4]]),
                              lambda u, v: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(share(v, dt[2], dt[3]))), share(u, dt[0], dt[1]) * math.sin(math.radians(share(v, dt[2], dt[3]))), dt[5]]))
            tfi.calc_grid(np.linspace(0, 1.0, dt[6]),
                          np.linspace(0, 1.0, dt[7]),
                          np.linspace(0, 1.0, dt[8]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_3d_sect-{}.xyz'.format(k))


"""
Implementation of the grid smoothing tools using elliptic PDE.

Note:
1. (i,j,k) is corresponding to (x,y,z),(u,v.w),(xi, eta, zeta),(x1,x2,x3),(I,J,K)...
2. (N,M) notation is not suggested to avoid confuse on column-and-row.
3. By default, delta_xi = delta_eta = delta_zeta = 1, and is neglected in code.
4. All the partial derivatives in elliptic PDE are calculated with the central scheme.
"""


class EllipticGrid2D(object):
    di = [0, 1, 0, -1, 0, 1, -1, -1, 1]
    dj = [0, 0, 1, 0, -1, 1, 1, -1, -1]

    def __init__(self, grid):
        """
        2D curvilinear grid based on the elliptic PDE.
        :param grid: Initial grid. The subscript iterate through (Dim1, Dim2, Dim3), each element contains (X, Y, Z).
        """

        '''Pre-check'''
        assert len(grid.shape) == 3
        assert grid.shape[-1] in (2, 3)

        '''Shape constants'''
        ii, jj, dim = grid.shape
        assert (ii - 2) * (jj - 2) != 0

        '''Grid'''
        self.r = np.copy(grid[:, :, :2])

        '''Partial derivatives'''
        self.r1 = np.zeros((ii, jj, 2))
        self.r2 = np.zeros((ii, jj, 2))
        self.r11 = np.zeros((ii, jj, 2))
        self.r22 = np.zeros((ii, jj, 2))
        self.r12 = np.zeros((ii, jj, 2))

        '''Source term'''
        self.pq = np.zeros((ii, jj, 2))

        '''Coefficients'''
        self.a = np.zeros((ii, jj))  # alpha
        self.b = np.zeros((ii, jj))  # beta
        self.g = np.zeros((ii, jj))  # gamma
        self.j2 = np.zeros((ii, jj))  # square of det(jacobi)

    @property
    def i_num(self):
        return self.r.shape[0]

    @property
    def j_num(self):
        return self.r.shape[1]

    def is_special(self, i, j):
        return i == 0 or i == self.i_num - 1 or j == 0 or j == self.j_num - 1

    def internal_pnt_idx(self, i, j):
        return (i - 1) + (j - 1) * (self.i_num - 2)

    @property
    def grid(self):
        return self.r

    def r_xi(self, i, j):
        return 0.5 * (self.r[i + 1][j] - self.r[i - 1][j])

    def r_eta(self, i, j):
        return 0.5 * (self.r[i][j + 1] - self.r[i][j - 1])

    def r_xi2(self, i, j):
        return self.r[i + 1][j] - 2 * self.r[i][j] + self.r[i - 1][j]

    def r_eta2(self, i, j):
        return self.r[i][j + 1] - 2 * self.r[i][j] + self.r[i][j - 1]

    def r_xi_eta(self, i, j):
        return 0.25 * (self.r[i + 1][j + 1] + self.r[i - 1][j - 1] - self.r[i + 1][j - 1] - self.r[i - 1][j + 1])

    def alpha(self, i, j):
        return vector_square(self.r_eta(i, j))

    def beta(self, i, j):
        return np.dot(self.r_eta(i, j), self.r_xi(i, j))

    def gamma(self, i, j):
        return vector_square(self.r_xi(i, j))

    def jacobi(self, i, j):
        return np.linalg.det(np.matrix([self.r_xi(i, j), self.r_eta(i, j)])) ** 2

    @abstractmethod
    def calc_all_param(self):
        pass

    @abstractmethod
    def calc_eqn_param(self, i, j):
        pass

    def smooth(self):
        """
        Smooth the grid with Picard iteration.
        :return: None.
        """

        var_num = (self.i_num - 2) * (self.j_num - 2)
        rhs = np.empty((var_num, 2))

        '''Solve the grid iteratively'''
        iteration_cnt = 0
        residual = sys.float_info.max
        while not math.isclose(residual, 0, abs_tol=1e-5):
            '''Calculate all coefficients and derivatives'''
            self.calc_all_param()

            '''Build Ax = b'''
            rhs.fill(0.0)
            eqn_idx = 0
            row = []
            col = []
            val = []
            for j in range(1, self.j_num - 1):
                for i in range(1, self.i_num - 1):
                    ca = self.calc_eqn_param(i, j)  # surrounding coefficients
                    for t in range(9):
                        ii = i + EllipticGrid2D.di[t]
                        jj = j + EllipticGrid2D.dj[t]
                        if self.is_special(ii, jj):
                            rhs[eqn_idx] -= ca[t] * self.r[ii][jj]
                        else:
                            row.append(eqn_idx)
                            col.append(self.internal_pnt_idx(ii, jj))
                            val.append(ca[t])
                    eqn_idx += 1

            '''Construct the sparse coefficient matrix'''
            scm = sparse.coo_matrix((val, (row, col)), shape=(var_num, var_num), dtype=float).tocsr()

            '''Solve the grid'''
            u = np.copy([dsolve.spsolve(scm, b) for b in rhs.transpose()]).transpose()

            '''Update and calculate residual'''
            eqn_idx = 0
            residual = sys.float_info.min
            for j in range(1, self.j_num - 1):
                for i in range(1, self.i_num - 1):
                    cur_residual = norm(u[eqn_idx] - self.r[i][j], np.inf)
                    if cur_residual > residual:
                        residual = cur_residual
                    self.r[i][j] = u[eqn_idx]
                    eqn_idx += 1

            iteration_cnt += 1


class Laplace2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        Smooth the grid with Laplace PDE.
        :param grid: The initial grid.
        """

        super(Laplace2D, self).__init__(grid)

    def calc_all_param(self):
        for i in range(1, self.i_num - 1):
            for j in range(1, self.j_num - 1):
                self.r1[i][j] = self.r_xi(i, j)
                self.r2[i][j] = self.r_eta(i, j)
                self.a[i][j] = vector_square(self.r2[i][j])
                self.b[i][j] = np.dot(self.r1[i][j], self.r2[i][j])
                self.g[i][j] = vector_square(self.r1[i][j])

    def calc_eqn_param(self, i, j):
        ans = np.empty(9, float)
        ans[0] = -2 * (self.a[i][j] + self.g[i][j])
        ans[1] = ans[3] = self.a[i][j]
        ans[2] = ans[4] = self.g[i][j]
        ans[6] = ans[8] = self.b[i][j] / 2
        ans[5] = ans[7] = -ans[6]
        return ans


class ThomasMiddlecoff2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        Smooth the grid with the Thomas-Middlecoff method.
        :param grid: Initial grid.
        """

        super(ThomasMiddlecoff2D, self).__init__(grid)

        self.phi = np.zeros((self.i_num, self.j_num))
        self.psi = np.zeros((self.i_num, self.j_num))

    def calc_all_param(self):
        for i in range(1, self.i_num - 1):
            for j in range(1, self.j_num - 1):
                self.r1[i][j] = self.r_xi(i, j)
                self.r2[i][j] = self.r_eta(i, j)
                self.r11[i][j] = self.r_xi2(i, j)
                self.r22[i][j] = self.r_eta2(i, j)
                self.r12[i][j] = self.r_xi_eta(i, j)
                self.a[i][j] = vector_square(self.r2[i][j])
                self.b[i][j] = np.dot(self.r1[i][j], self.r2[i][j])
                self.g[i][j] = vector_square(self.r1[i][j])
        for j in (0, self.j_num - 1):
            for i in range(1, self.i_num - 1):
                self.r1[i][j] = self.r_xi(i, j)
                self.r11[i][j] = self.r_xi2(i, j)
                self.phi[i][j] = - np.dot(self.r1[i][j], self.r11[i][j]) / vector_square(self.r1[i][j])
        for i in (0, self.i_num - 1):
            for j in range(1, self.j_num - 1):
                self.r2[i][j] = self.r_eta(i, j)
                self.r22[i][j] = self.r_eta2(i, j)
                self.psi[i][j] = -np.dot(self.r2[i][j], self.r22[i][j]) / vector_square(self.r2[i][j])
        for i in range(1, self.i_num - 1):
            dist = np.linspace(self.phi[i][0], self.phi[i][self.j_num - 1], self.j_num)
            for j in range(1, self.j_num - 1):
                self.phi[i][j] = dist[j]
        for j in range(1, self.j_num - 1):
            dist = np.linspace(self.psi[0][j], self.psi[self.i_num - 1][j], self.i_num)
            for i in range(1, self.i_num - 1):
                self.psi[i][j] = dist[i]

    def calc_eqn_param(self, i, j):
        ans = np.empty(9, float)
        ans[0] = -2.0 * (self.a[i][j] + self.g[i][j])
        ans[1] = self.a[i][j] * (1 + self.phi[i][j] / 2)
        ans[2] = self.g[i][j] * (1 + self.psi[i][j] / 2)
        ans[3] = self.a[i][j] * (1 - self.phi[i][j] / 2)
        ans[4] = self.g[i][j] * (1 - self.psi[i][j] / 2)
        ans[6] = ans[8] = self.b[i][j] / 2
        ans[5] = ans[7] = -ans[6]
        return ans


class EllipticGrid3D(object):
    di = [0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1]
    dj = [0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1, 0, 0, 0, 0]
    dk = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1]

    def __init__(self, grid):
        """
        Base class for 3D elliptic-PDE based grid smoother.
        :param grid: Initial grid.
        """

        '''Pre-check'''
        assert len(grid.shape) == 4
        assert grid.shape[-1] == 3

        '''Shape constants'''
        ii, jj, kk, dim = grid.shape
        assert (ii - 2) * (jj - 2) * (kk - 2) != 0

        '''Grid'''
        self.r = np.copy(grid)

        '''Derivatives'''
        self.r1 = np.zeros_like(self.r)
        self.r2 = np.zeros_like(self.r)
        self.r3 = np.zeros_like(self.r)
        self.r11 = np.zeros_like(self.r)
        self.r22 = np.zeros_like(self.r)
        self.r33 = np.zeros_like(self.r)
        self.r12 = np.zeros_like(self.r)
        self.r23 = np.zeros_like(self.r)
        self.r31 = np.zeros_like(self.r)

        '''Source term'''
        self.pqr = np.zeros((ii, jj, kk, 3))

        '''Equation coefficients'''
        self.a1 = np.zeros((ii, jj, kk))  # alpha1
        self.a2 = np.zeros((ii, jj, kk))  # alpha2
        self.a3 = np.zeros((ii, jj, kk))  # alpha3
        self.b12 = np.zeros((ii, jj, kk))  # beta12
        self.b23 = np.zeros((ii, jj, kk))  # beta23
        self.b31 = np.zeros((ii, jj, kk))  # beta23
        self.j2 = np.zeros((ii, jj, kk))  # Jacobi

    @property
    def i_num(self):
        return self.r.shape[0]

    @property
    def j_num(self):
        return self.r.shape[1]

    @property
    def k_num(self):
        return self.r.shape[2]

    def is_special(self, i, j, k):
        return i == 0 or i == self.i_num - 1 or j == 0 or j == self.j_num - 1 or k == 0 or k == self.k_num - 1

    def internal_pnt_idx(self, i, j, k):
        return (k - 1) * (self.j_num - 2) * (self.i_num - 2) + (j - 1) * (self.i_num - 2) + (i - 1)

    @property
    def grid(self):
        return self.r

    def r_xi(self, i, j, k):
        return 0.5 * (self.r[i + 1][j][k] - self.r[i - 1][j][k])

    def r_eta(self, i, j, k):
        return 0.5 * (self.r[i][j + 1][k] - self.r[i][j - 1][k])

    def r_zeta(self, i, j, k):
        return 0.5 * (self.r[i][j][k + 1] - self.r[i][j][k - 1])

    def r_xi2(self, i, j, k):
        return self.r[i - 1][j][k] - 2 * self.r[i][j][k] + self.r[i + 1][j][k]

    def r_eta2(self, i, j, k):
        return self.r[i][j - 1][k] - 2 * self.r[i][j][k] + self.r[i][j + 1][k]

    def r_zeta2(self, i, j, k):
        return self.r[i][j][k - 1] - 2 * self.r[i][j][k] + self.r[i][j][k + 1]

    def r_xi_eta(self, i, j, k):
        return 0.25 * (self.r[i + 1][j + 1][k] - self.r[i + 1][j - 1][k] - self.r[i - 1][j + 1][k] + self.r[i - 1][j - 1][k])

    def r_eta_zeta(self, i, j, k):
        return 0.25 * (self.r[i][j + 1][k + 1] - self.r[i][j - 1][k + 1] - self.r[i][j + 1][k - 1] + self.r[i][j - 1][k - 1])

    def r_zeta_xi(self, i, j, k):
        return 0.25 * (self.r[i + 1][j][k + 1] - self.r[i - 1][j][k + 1] - self.r[i + 1][j][k - 1] + self.r[i - 1][j][k - 1])

    def alpha1(self, i, j, k):
        r2 = self.r_eta(i, j, k)
        r3 = self.r_zeta(i, j, k)
        return vector_square(r2) * vector_square(r3) - np.dot(r2, r3) ** 2

    def alpha2(self, i, j, k):
        r3 = self.r_zeta(i, j, k)
        r1 = self.r_xi(i, j, k)
        return vector_square(r3) * vector_square(r1) - np.dot(r3, r1) ** 2

    def alpha3(self, i, j, k):
        r1 = self.r_xi(i, j, k)
        r2 = self.r_eta(i, j, k)
        return vector_square(r1) * vector_square(r2) - np.dot(r1, r2) ** 2

    def beta12(self, i, j, k):
        r1 = self.r_xi(i, j, k)
        r2 = self.r_eta(i, j, k)
        r3 = self.r_zeta(i, j, k)
        return np.dot(r1, r3) * np.dot(r2, r3) - np.dot(r1, r2) * norm(r3)

    def beta23(self, i, j, k):
        r1 = self.r_xi(i, j, k)
        r2 = self.r_eta(i, j, k)
        r3 = self.r_zeta(i, j, k)
        return np.dot(r2, r1) * np.dot(r3, r1) - np.dot(r2, r3) * norm(r1)

    def beta31(self, i, j, k):
        r1 = self.r_xi(i, j, k)
        r2 = self.r_eta(i, j, k)
        r3 = self.r_zeta(i, j, k)
        return np.dot(r3, r2) * np.dot(r1, r2) - np.dot(r3, r1) * norm(r2)

    def jacobi(self, i, j, k):
        return np.linalg.det(np.matrix([self.r_xi(i, j, k), self.r_eta(i, j, k), self.r_zeta(i, j, k)])) ** 2

    @abstractmethod
    def calc_all_param(self):
        pass

    @abstractmethod
    def calc_eqn_param(self, i, j, k):
        pass

    def smooth(self):
        """
        Smooth the grid using Picard Iteration.
        :return: None.
        """

        var_num = (self.i_num - 2) * (self.j_num - 2) * (self.k_num - 2)
        rhs = np.zeros((var_num, 3))

        '''Solve the grid iteratively'''
        iteration_cnt = 0
        residual = sys.float_info.max
        while not math.isclose(residual, 0, abs_tol=1e-5):
            '''Calculate all coefficients'''
            self.calc_all_param()

            '''Build Ax=b'''
            eqn_idx = 0
            rhs.fill(0.0)
            row = []
            col = []
            val = []
            for k in range(1, self.k_num - 1):
                for j in range(1, self.j_num - 1):
                    for i in range(1, self.i_num - 1):
                        '''Calculate stencil coefficients'''
                        ca = self.calc_eqn_param(i, j, k)

                        '''Construct the equation'''
                        for t in range(19):
                            ii = i + EllipticGrid3D.di[t]
                            jj = j + EllipticGrid3D.dj[t]
                            kk = k + EllipticGrid3D.dk[t]
                            if self.is_special(ii, jj, kk):
                                rhs[eqn_idx] -= ca[t] * self.r[ii][jj][kk]
                            else:
                                row.append(eqn_idx)
                                col.append(self.internal_pnt_idx(ii, jj, kk))
                                val.append(ca[t])

                        eqn_idx += 1

            '''Construct the sparse coefficient matrix'''
            scm = sparse.coo_matrix((val, (row, col)), shape=(var_num, var_num), dtype=float).tocsr()

            '''Solve the grid'''
            u = np.copy([dsolve.spsolve(scm, b) for b in rhs.transpose()]).transpose()

            '''Update and calculate residual'''
            eqn_idx = 0
            residual = sys.float_info.min
            for k in range(1, self.k_num - 1):
                for j in range(1, self.j_num - 1):
                    for i in range(1, self.i_num - 1):
                        cur_residual = norm(u[eqn_idx] - self.r[i][j][k], np.inf)
                        if cur_residual > residual:
                            residual = cur_residual
                        self.r[i][j][k] = u[eqn_idx]
                        eqn_idx += 1

            iteration_cnt += 1
            print(residual)


class Laplace3D(EllipticGrid3D):
    def __init__(self, grid):
        super(Laplace3D, self).__init__(grid)

    def calc_all_param(self):
        for i in range(1, self.i_num - 1):
            for j in range(1, self.j_num - 1):
                for k in range(1, self.k_num - 1):
                    self.r1[i][j][k] = self.r_xi(i, j, k)
                    self.r2[i][j][k] = self.r_eta(i, j, k)
                    self.r3[i][j][k] = self.r_zeta(i, j, k)
                    self.r11[i][j][k] = self.r_xi2(i, j, k)
                    self.r22[i][j][k] = self.r_eta2(i, j, k)
                    self.r33[i][j][k] = self.r_zeta2(i, j, k)
                    self.r12[i][j][k] = self.r_xi_eta(i, j, k)
                    self.r23[i][j][k] = self.r_eta_zeta(i, j, k)
                    self.r31[i][j][k] = self.r_zeta_xi(i, j, k)
                    self.a1[i][j][k] = vector_square(self.r2[i][j][k]) * vector_square(self.r3[i][j][k]) - np.dot(self.r2[i][j][k], self.r3[i][j][k]) ** 2
                    self.a2[i][j][k] = vector_square(self.r3[i][j][k]) * vector_square(self.r1[i][j][k]) - np.dot(self.r3[i][j][k], self.r1[i][j][k]) ** 2
                    self.a3[i][j][k] = vector_square(self.r1[i][j][k]) * vector_square(self.r2[i][j][k]) - np.dot(self.r1[i][j][k], self.r2[i][j][k]) ** 2
                    self.b12[i][j][k] = np.dot(self.r1[i][j][k], self.r3[i][j][k]) * np.dot(self.r2[i][j][k], self.r3[i][j][k]) - np.dot(self.r1[i][j][k], self.r2[i][j][k]) * norm(self.r3[i][j][k])
                    self.b23[i][j][k] = np.dot(self.r2[i][j][k], self.r1[i][j][k]) * np.dot(self.r3[i][j][k], self.r1[i][j][k]) - np.dot(self.r2[i][j][k], self.r3[i][j][k]) * norm(self.r1[i][j][k])
                    self.b31[i][j][k] = np.dot(self.r3[i][j][k], self.r2[i][j][k]) * np.dot(self.r1[i][j][k], self.r2[i][j][k]) - np.dot(self.r3[i][j][k], self.r1[i][j][k]) * norm(self.r2[i][j][k])

    def calc_eqn_param(self, i, j, k):
        ans = np.empty(19, float)
        ans[0] = -2 * (self.a1[i][j][k] + self.a2[i][j][k] + self.a3[i][j][k])
        ans[1] = ans[2] = self.a1[i][j][k]
        ans[3] = ans[4] = self.a2[i][j][k]
        ans[5] = ans[6] = self.a3[i][j][k]
        ans[7] = ans[8] = 0.5 * self.b12[i][j][k]
        ans[9] = ans[10] = -ans[8]
        ans[11] = ans[12] = 0.5 * self.b23[i][j][k]
        ans[13] = ans[14] = -ans[12]
        ans[15] = ans[16] = 0.5 * self.b31[i][j][k]
        ans[17] = ans[18] = -ans[16]
        return ans


class ThomasMiddlecoff3D(EllipticGrid3D):
    def __init__(self, grid):
        super(ThomasMiddlecoff3D, self).__init__(grid)

        self.phi = np.zeros((self.i_num, self.j_num, self.k_num))
        self.psi = np.zeros((self.i_num, self.j_num, self.k_num))
        self.omega = np.zeros((self.i_num, self.j_num, self.k_num))

    def calc_all_param(self):
        pass

    def calc_eqn_param(self, i, j, k):
        pass


class EllipticTestCase(unittest.TestCase):
    def test_3d_laplace(self):
        # Delta, R1, R2, U, V, W
        data = [(-10, 4, 25, 16, 41, 21),
                (-1, 2, 5, 16, 33, 16)]

        tfi_grid = []
        laplace_grid = []
        for dt in data:
            tfi = LinearTFI2D(lambda u: np.array([(dt[0] + dt[1]) * (1 - u) + dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[1] * math.cos(math.pi * v) + dt[0], dt[1] * math.sin(math.pi * v), 0]),
                              lambda u: np.array([(dt[0] - dt[1]) * (1 - u) - dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[2] * math.cos(math.pi * v), dt[2] * math.sin(math.pi * v), 0]))
            tfi.calc_grid(uniform(dt[3]), uniform(dt[4]))
            extruded_grid = np.empty((dt[3], dt[4], dt[5], 3))
            for i in range(dt[3]):
                for j in range(dt[4]):
                    for k in range(dt[5]):
                        extruded_grid[i][j][k] = tfi.grid[i][j]
                        extruded_grid[i][j][k][2] = k
            laplace_blk = Laplace3D(extruded_grid)
            laplace_blk.smooth()
            tfi_grid.append(extruded_grid)
            laplace_grid.append(laplace_blk.grid)

        for i in range(len(data)):
            msh = Plot3D()
            msh.add(Plot3DBlock.construct_from_array(tfi_grid[i]))
            msh.save('test_3d_eccentric-{}_tfi.xyz'.format(i))
            msh.clear()
            msh.add(Plot3DBlock.construct_from_array(laplace_grid[i]))
            msh.save('test_3d_eccentric-{}_laplace.xyz'.format(i))
        self.assertTrue(True)


"""
Implementation of the Plot3D standard.

Note:
All the units are SI by default.
"""


class Plot3DBlock(object):
    def __init__(self, pts):
        """
        Single Block structured grid.
        IBLANK Values:
             0 : Outside of the computational area.
             1 : Normal points.
             2 : Wall boundary points.
            -n : Adjacent to the n-th block.
        :param pts: All the coordinates, the index traverse through(I, J, K) from left to right, each coordinate is consist of (X, Y, Z, IBLANK).
        """

        self.data = np.copy(pts)

    @property
    def dim(self):
        """
        Dimension of the grid.
        """

        return self.data.shape[:3]

    def __repr__(self):
        ii, jj, kk = self.dim
        return "{} x {}".format(ii, jj) if kk == 1 else "{} x {} x {}".format(ii, jj, kk)

    @property
    def pnt_num(self):
        """
        :return Number of points.
        :rtype int
        """

        ii, jj, kk = self.dim
        return ii * jj * kk

    @property
    def cell_num(self):
        """
        :return Number of cells.
        :rtype int
        """

        ii, jj, kk = self.dim
        return (ii - 1) * (jj - 1) if kk == 1 else (ii - 1) * (jj - 1) * (kk - 1)

    @property
    def face_num(self):
        """
        :return Number of faces.
        :rtype int
        """

        ii, jj, kk = self.dim
        return 3 * ii * jj * kk - (ii * jj + jj * kk + kk * ii)

    @property
    def all_pnt(self):
        t = 0
        ii, jj, kk = self.dim
        ret = np.empty((self.pnt_num, 3), float)
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    ret[t] = self.data[i][j][k][:3]
                    t += 1

        return ret

    def write(self, f_out, with_iblank):
        """
        Output the grid into a stream.
        :param f_out: The output stream.
        :param with_iblank: Indicate if the IBLANK info is included.
        :type with_iblank: bool
        :return: None
        """

        ii, jj, kk = self.dim
        for d in range(3):
            for k in range(kk):
                for j in range(jj):
                    for i in range(ii):
                        f_out.write("{}{}".format('\n' if i == 0 else ' ', self.data[i][j][k][d]))

        if with_iblank:
            for k in range(kk):
                for j in range(jj):
                    for i in range(ii):
                        f_out.write("{}{}".format('\n' if i == 0 else ' ', int(self.data[i][j][k][-1])))

    def set_iblank(self, i, j, k, t):
        """
        Set the IBLANK value on certain point.
        :param i: Index in X direction.
        :type i: int
        :param j: Index in Y direction.
        :type j: int
        :param k: Index in Z direction.
        :type k: int
        :param t: Target IBLANK Value.
        :type t: int
        :return: None.
        """

        self.data[i][j][k][-1] = t

    def set_area_iblank(self, rg0, rg1, rg2, t):
        """
        设置区域内网格点的IBLANK信息
        :param rg0: X方向范围
        :param rg1: Y方向范围
        :param rg2: Z方向范围
        :param t: IBLANK Value
        :type t: int
        :return: None
        """

        for i in rg0:
            for j in rg1:
                for k in rg2:
                    self.set_iblank(i, j, k, t)

    def set_boundary_iblank(self, t):
        """
        设置边界上网格点的IBLANK信息
        :param t: IBLANK Value
        :type t: int
        :return: None
        """

        ii, jj, kk = self.dim
        if kk == 1:
            for i in range(ii):
                for j in range(jj):
                    if i in (0, ii - 1) or j in (0, jj - 1):
                        self.set_iblank(i, j, 0, t)
        else:
            for i in range(ii):
                for j in range(jj):
                    for k in range(kk):
                        if i in (0, ii - 1) or j in (0, jj - 1) or k in (0, kk - 1):
                            self.set_iblank(i, j, k, t)

    @classmethod
    def construct_from_array(cls, pts):
        """
        Construct the Plot3D Block from grid array.
        :param pts: Input grid points, the index traverse (I, J)/(I, J, K) from left to right,
                    point should be 3D, each element is consist of (X, Y, Z).
        :return: Single-Block grid in Plot3D notation.
        :rtype: Plot3DBlock
        """

        if len(pts.shape) == 3:
            ni, nj, nd = pts.shape
            nk = 1
        elif len(pts.shape) == 4:
            ni, nj, nk, nd = pts.shape
        else:
            raise AssertionError("Invalid input grid array.")

        p3d = np.zeros((ni, nj, nk, 4))
        if len(pts.shape) == 3:
            for i in range(ni):
                for j in range(nj):
                    for d in range(nd):
                        p3d[i][j][0][d] = pts[i][j][d]
                    p3d[i][j][0][-1] = 1
        else:
            for i in range(ni):
                for j in range(nj):
                    for k in range(nk):
                        for d in range(nd):
                            p3d[i][j][k][d] = pts[i][j][k][d]
                        p3d[i][j][k][-1] = 1

        ret = cls(p3d)
        ret.set_boundary_iblank(2)  # Close the boundary by default.
        return ret


class Plot3D(object):
    def __init__(self):
        self.blk_list = []

    def __repr__(self):
        ret = "Multi-Block structured grid in Plot3D format with {} block(s)\n".format(self.size)
        for i in range(self.size):
            ret += "{}: {}\n".format(i, repr(self.blk_list[i]))

        return ret

    @property
    def size(self):
        return len(self.blk_list)

    def add(self, blk):
        """
        Append new block.
        :param blk: single block.
        :type blk: Plot3DBlock
        :return: None
        """

        self.blk_list.append(blk)

    def clear(self):
        self.blk_list.clear()

    def save(self, fn, with_iblank=False):
        """
        Output grid into file.
        :param fn: Filename.
        :type fn: str
        :param with_iblank: Indicate if the IBLANK info is included.
        :type with_iblank: bool
        :return: None
        """

        f_out = open(fn, 'w')
        f_out.write("{}".format(self.size))

        for blk in self.blk_list:
            ii, jj, kk = blk.dim
            f_out.write("\n{} {} {}".format(ii, jj, kk))

        for k, blk in enumerate(self.blk_list):
            blk.write(f_out, with_iblank)

        f_out.close()


class Plot3DTestCase(unittest.TestCase):
    def test_single(self):
        # x_min, x_max, y_min, y_max, z_min, z_max, nu, nv, nw
        rect_param = [(0, 100, 0, 60, 0, 40, 61, 16, 21),
                      (0, 100, 0, 60, 0, 40, 61, 16, 1)]
        # r_min, r_max, theta_min, theta_max, h_min, h_max, nu, nv, nw
        sect_param = [(50, 100, 60, 320, 0, 30, 61, 16, 21),
                      (50, 100, 60, 320, 0, 30, 61, 16, 1)]
        ans = []

        for p in rect_param:
            x_min, x_max, y_min, y_max, z_min, z_max, nu, nv, nw = p
            pts = np.zeros((nu, nv, nw, 3))
            u_list = np.linspace(x_min, x_max, nu)
            v_list = np.linspace(y_min, y_max, nv)
            w_list = np.linspace(z_min, z_max, nw)
            for i in range(nu):
                for j in range(nv):
                    for k in range(nw):
                        pts[i][j][k][0] = u_list[i]
                        pts[i][j][k][1] = v_list[j]
                        pts[i][j][k][2] = w_list[k]
            ans.append(pts)

        for p in sect_param:
            r_min, r_max, theta_min, theta_max, h_min, h_max, nu, nv, nw = p
            pts = np.zeros((nu, nv, nw, 3))
            u_list = np.linspace(r_min, r_max, nu)
            v_list = np.linspace(theta_min, theta_max, nv)
            w_list = np.linspace(h_min, h_max, nw)
            for i in range(nu):
                for j in range(nv):
                    for k in range(nw):
                        ct = math.radians(v_list[j])
                        pts[i][j][k] = np.array([u_list[i] * math.cos(ct), u_list[i] * math.sin(ct), w_list[k]])
            ans.append(pts)

        self.assertTrue(len(ans) == len(rect_param) + len(sect_param))

        grid = Plot3D()
        for t in range(len(ans)):
            grid.clear()
            blk = Plot3DBlock.construct_from_array(ans[t])
            grid.add(blk)
            fn = "'test_plot3d_single-{}.xyz".format(t)
            grid.save(fn)

    def test_multi(self):
        # x_min, x_max, y_min, y_max, z_min, z_max, nu, nv, nw
        rect_param = [(0, 100, 0, 60, 0, 40, 61, 16, 21),
                      (120, 200, 75, 100, 50, 80, 61, 16, 21)]
        ans = []

        for p in rect_param:
            x_min, x_max, y_min, y_max, z_min, z_max, nu, nv, nw = p
            pts = np.zeros((nu, nv, nw, 3))
            u_list = np.linspace(x_min, x_max, nu)
            v_list = np.linspace(y_min, y_max, nv)
            w_list = np.linspace(z_min, z_max, nw)
            for i in range(nu):
                for j in range(nv):
                    for k in range(nw):
                        pts[i][j][k][0] = u_list[i]
                        pts[i][j][k][1] = v_list[j]
                        pts[i][j][k][2] = w_list[k]
            ans.append(pts)

        self.assertTrue(len(ans) == len(rect_param))

        grid = Plot3D()
        for t in range(len(ans)):
            blk = Plot3DBlock.construct_from_array(ans[t])
            grid.add(blk)
        grid.save('test_plot3d_multi.xyz')


"""
Implementation of NASA's Neutral Map File representation.

Ref:
(1) https://geolab.larc.nasa.gov/Volume/Doc/nmf.htm
(2) https://github.com/mrklein/p3d2gmsh/blob/master/p3d2gmsh.py
"""


class NMFEntry(object):
    NMF_LITERAL = '#{:^19}{:^8}{:^2}{:^8}{:^8}{:^8}{:^8}{:^8}{:^2}{:^8}{:^8}{:^8}{:^8}{:^6}'.format('Type', 'B1', 'F1', 'S1', 'E1', 'S2', 'E2', 'B2', 'F2', 'S1', 'E1', 'S2', 'E2', 'Swap')

    def __init__(self, tp, b1, f1, rg1, b2, f2, rg2, swp):
        """
        Entry used in Neutral Map File to indicate the topological features of the mesh.
        :param tp:The type of feature (topological or boundary condition) to be defined and positioned within the mesh.
        :type tp: str
        :param b1: The number of the first block(Starting from 1).
        :type b1: int
        :param f1: The face number for the first block(From 1 to 6).
        :type f1: int
        :param rg1: Range descriptions of the first block.
        :param b2:The number of the second block(Starting from 1).
        :type b2: int
        :param f2: The face number for the second block(From 1 to 6).
        :type f2: int
        :param rg2: Range descriptions of the second block.
        :param swp: Orientation flag(Specified only for Type==ONE_TO_ONE).
                    False if the primary directions of the two identified faces are aligned(though perhaps in opposite directions) and True otherwise.
        :type swp: bool
        """

        assert b1 > 0 and 1 <= f1 <= 6 and rg1.shape == (2, 2)
        assert b2 > 0 and 1 <= f2 <= 6 and rg2.shape == (2, 2)

        self.Type = tp
        self.B1 = b1
        self.F1 = f1
        self.B1PriStart = rg1[0][0]  # The starting index in the primary coordinate direction for the face in the 1st block.
        self.B1PriEnd = rg1[0][1]  # The ending index in the primary coordinate direction for the face in the 1st block.
        self.B1SecStart = rg1[1][0]  # The starting index in the secondary coordinate direction for the face in the 1st block.
        self.B1SecEnd = rg1[1][1]  # The ending index in the secondary coordinate direction for the face in the 1st block.
        self.B2 = b2
        self.F2 = f2
        self.B2PriStart = rg2[0][0]  # The starting index in the primary coordinate direction for the face in the 2nd block.
        self.B2PriEnd = rg2[0][1]  # The ending index in the primary coordinate direction for the face in the 2nd block.
        self.B2SecStart = rg2[1][0]  # The starting index in the secondary coordinate direction for the face in the 2nd block.
        self.B2SecEnd = rg2[1][1]  # The ending index in the secondary coordinate direction for the face in the 2nd block.
        self.Swap = swp

    def write(self, f_out):
        ret = '{:<20}'.format(self.Type)
        ret += '{:>8}'.format(self.B1)
        ret += '{:>2}'.format(self.F1)
        ret += '{:>8}'.format(self.B1PriStart)
        ret += '{:>8}'.format(self.B1PriEnd)
        ret += '{:>8}'.format(self.B1SecStart)
        ret += '{:>8}'.format(self.B1SecEnd)
        ret += '{:>8}'.format(self.B2)
        ret += '{:>2}'.format(self.F2)
        ret += '{:>8}'.format(self.B2PriStart)
        ret += '{:>8}'.format(self.B2PriEnd)
        ret += '{:>8}'.format(self.B2SecStart)
        ret += '{:>8}'.format(self.B2SecEnd)
        ret += '{:>6}'.format('TRUE' if self.Swap else 'FALSE')
        f_out.write(ret)


def blk_node_num(shape):
    """
    Calculate the num of nodes in a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of nodes in the block.
    :rtype: int
    """

    assert len(shape) in (3, 4)
    if len(shape) == 3:
        u, v, _ = shape
        return u * v
    else:
        u, v, w, _ = shape
        return u * v * w


def blk_face_num(shape):
    """
    Calculate the num of faces in a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of faces in the block.
    :rtype: int
    """

    assert len(shape) in (3, 4)
    if len(shape) == 3:
        u, v, _ = shape
        return u * (v - 1) + v * (u - 1)
    else:
        u, v, w, _ = shape
        return u * (v - 1) * (w - 1) + v * (w - 1) * (u - 1) + w * (u - 1) * (v - 1)


def blk_cell_num(shape):
    """
    Calculate the num of cells in a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of cells in the block.
    :rtype: int
    """

    assert len(shape) in (3, 4)
    if len(shape) == 3:
        u, v, _ = shape
        return (u - 1) * (v - 1)
    else:
        u, v, w, _ = shape
        return (u - 1) * (v - 1) * (w - 1)


def blk_internal_node_num(shape):
    """
    Calculate the num of nodes inside a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of cells inside the block.
    :rtype: int
    """

    assert len(shape) in (3, 4)
    if len(shape) == 3:
        u, v, _ = shape
        return (u - 2) * (v - 2)
    else:
        u, v, w, _ = shape
        return (u - 2) * (v - 2) * (w - 2)


def blk_internal_face_num(shape):
    """
    Calculate the num of faces inside a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of faces inside the block.
    :rtype: int
    """

    assert len(shape) in (3, 4)
    if len(shape) == 3:
        u, v, _ = shape
        return (u - 2) * (v - 1) + (u - 1) * (v - 2)
    else:
        u, v, w, _ = shape
        return (u - 1) * (v - 1) * (w - 2) + (u - 1) * (v - 2) * (w - 1) + (u - 2) * (v - 1) * (w - 1)


def on_boundary(pnt, shape):
    """
    Judge if the point is on the boundary of the block.
    :param pnt: Coordinate of the point.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: If the point is on the boundary.
    :rtype: bool
    """

    assert len(pnt) == len(shape) - 1
    assert len(shape) in (3, 4)

    if len(shape) == 3:
        i, j = pnt
        u, v, _ = shape
        return i == 0 or i == u - 1 or j == 0 or j == v - 1
    else:
        i, j, k = pnt
        u, v, w, _ = shape
        return i == 0 or i == u - 1 or j == 0 or j == v - 1 or k == 0 or j == _w - 1


class NeutralMapFile(object):
    def __init__(self, str_grid):
        self.blk = str_grid
        self.desc = []

        self.shell_pnt_idx = []
        self.shell_pnt_num = np.empty(blk_num, int)
        for k in range(blk_num):
            pn = boundary_node_num(blk_shape[k])
            shell_pnt_num[k] = pn
            shell_pnt_idx.append(np.zeros(pn, int))

    @property
    def grid(self):
        return self.blk

    @property
    def dim(self):
        return len(self.blk[0].shape) - 1

    @property
    def blk_num(self):
        return len(self.blk)

    @property
    def cell_num(self):
        return sum([blk_cell_num(self.blk[i].shape) for i in range(self.blk_num)])

    def add(self, entry):
        """
        Add topological or boundary info of the grid.
        :param entry: Description of the topological info.
        :type entry: NMFEntry
        :return: None.
        """

        self.desc.append(entry)

    def calc_topology(self):
        pass

    def save(self, fn):
        """
        Save the boundary mapping info into file.
        :param fn: File name.
        :type fn: str
        :return: None.
        """

        f_out = open(fn, 'w')
        f_out.write('# =========================== Neutral Map File generated by the Sky software =================================')
        f_out.write('# ============================================================================================================')
        f_out.write('# Block#    IDIM    JDIM    KDIM')
        f_out.write('# ------------------------------------------------------------------------------------------------------------')
        f_out.write('{:>8}'.format(self.blk_num))
        for i in range(self.blk_num):
            u, v, w, _ = self.blk[i].shape
            f_out.write('{:>8}{:>8}{:>8}{:>8}'.format(i + 1, u, v, w))
        f_out.write('# ============================================================================================================')
        f_out.write(NMFEntry.NMF_LITERAL)
        f_out.write('# ------------------------------------------------------------------------------------------------------------')
        for i in range(len(self.desc)):
            self.desc[i].write(f_out)
        f_out.close()


class NMFTestCase(unittest.TestCase):
    def test_pnt_num(self):
        pass


"""
Implementation of the ANSYS Fluent MSH File Format.

Note:
Refer from the appendix of ANSYS Fluent 15.0 User Manual. 
"""


class XFSection(object):
    def __init__(self, idx):
        """
        Abstract class of sections in the ANSYS Fluent MSH File Format.
        :param idx: Category of the section.
        :type idx: int
        """

        self.index = idx
        self.formatted_content = None

    @abstractmethod
    def build_content(self):
        pass

    @property
    def content(self):
        return self.formatted_content

    @abstractmethod
    def write(self, out_stream):
        pass


class XFComment(XFSection):
    def __init__(self, msg=''):
        """
        Comment Section in the ANSYS Fluent MSH File Format.
        :param msg: The comment message.
        :type msg: str
        """

        super(XFComment, self).__init__(0)
        self.msg = msg

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.msg)

    def write(self, out_stream):
        out_stream.write("({} \"{}\")".format(self.index, self.msg))


class XFHeader(XFSection):
    def __init__(self, info=''):
        """
        Header Section in the ANSYS Fluent MSH File Format.
        It's used to identify the program that wrote the file.
        :param info: Header info.
        :type info: str
        """

        super(XFHeader, self).__init__(1)
        self.info = "Grid generated with Python " + platform.python_version() if info == '' else info

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.info)

    def write(self, out_stream):
        out_stream.write("({} \"{}\")".format(self.index, self.info))


class XFDimension(XFSection):
    def __init__(self, dim):
        """
        Dimension Section in the ANSYS Fluent MSH File Format.
        It's used to specified the dimension of the grid.
        :param dim: Dimension
        :type dim: int
        """

        assert dim in (2, 3)
        super(XFDimension, self).__init__(2)
        self.ND = dim

    def build_content(self):
        self.formatted_content = "({} {})".format(self.index, self.ND)

    def write(self, out_stream):
        out_stream.write("({} {})".format(self.index, self.ND))


@unique
class NodeType(Enum):
    Virtual, Any, Boundary = range(3)


class XFNode(XFSection):
    def __init__(self, zone, first, last, tp, dim, pts=None):
        """
        Node Section in the ANSYS Fluent MSH File Format.
        :param zone: 区域号
        :type zone: int
        :param first: 网格点起始序号(Starting from 1)
        :type first: int
        :param last:  网格点终止序号
        :type last: int
        :param tp: 网格点类型(TGrid usage)
        :type tp: NodeType
        :param dim: Indicate the dimensionality of the node data.
        :type dim: int
        :param pts: 网格点数组，(last-first+1)个元素
        """

        super(XFNode, self).__init__(10)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.node_type = tp.value
        self.ND = dim
        self.pts = None if pts is None else np.copy(pts)

    def build_content(self):
        if self.pts is None:
            self.formatted_content = "({} ({} {} {} {})".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.node_type)[2:])
        else:
            self.formatted_content = "({} ({} {} {} {} {})".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.node_type)[2:], self.ND)
            n, dim = self.pts.shape
            self.formatted_content += "("
            for i in range(n):
                for d in range(self.ND):
                    self.formatted_content += "{}{}".format('\n' if d == 0 else ' ', self.pts[i][d])
            self.formatted_content += ")"
        self.formatted_content += ")"

    @classmethod
    def declaration(cls, num):
        return cls(0, 1, num, NodeType.Virtual, 0)

    def write(self, out_stream):
        if self.pts is None:
            out_stream.write("({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.node_type)[2:]))
        else:
            out_stream.write("({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.node_type)[2:], self.ND))
            n, dim = self.pts.shape
            for i in range(n):
                out_stream.write("\n{}".format(self.pts[i][0]))
                for d in range(1, self.ND):
                    out_stream.write(" {}".format(self.pts[i][d]))
            out_stream.write("))")


class BCType(Enum):
    Inactive = 0
    Interior = 2
    Wall = 3
    PressureInlet = 4
    InletVent = 4
    IntakeFan = 4
    PressureOutlet = 5
    OutletVent = 5
    ExhaustFan = 5
    Symmetry = 7
    PeriodicShadow = 8
    PressureFarField = 9
    VelocityInlet = 10
    Periodic = 12
    Fan = 14
    PorousJump = 14
    Radiator = 14
    MassFlowInlet = 20
    Interface = 24
    Parent = 31  # hanging node
    Outflow = 36
    Axis = 37


@unique
class FaceType(Enum):
    Mixed = 0
    Linear = 2
    Triangular = 3
    Quadrilateral = 4
    Polygonal = 5


class XFFace(XFSection):
    def __init__(self, zone, first, last, bct, ft, face_info=None):
        """
        Face Section in the ANSYS Fluent MSH File Format.
        :param zone: 区域号
        :type zone: int
        :param first: 起始序号(Starting from 1)
        :type first: int
        :param last: 终止序号
        :type last: int
        :param bct: 边界属性
        :type bct: BCType
        :param ft: 形状类别
        :type ft: FaceType
        :param face_info: 边界信息，一维数组(last-first+1)个元素，每个元素中包含了邻接关系
        """

        super(XFFace, self).__init__(13)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.bc_type = bct.value
        self.face_type = ft.value
        self.face_info = None if face_info is None else np.copy(face_info)

    def build_content(self):
        if self.face_info is not None:
            self.formatted_content = "({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.bc_type)[2:], hex(self.face_type)[2:])
            for fc in self.face_info:
                for cfi in range(len(fc)):  # current face info
                    self.formatted_content += "{}{}".format('\n' if cfi == 0 else ' ', hex(fc[cfi])[2:])
            self.formatted_content += '))'
        else:
            self.formatted_content = "({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.face_type)[2:])

    @classmethod
    def declaration(cls, num):
        return cls(0, 1, num, BCType.Inactive, FaceType.Mixed)

    def write(self, out_stream):
        if self.face_info is not None:
            out_stream.write("({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.bc_type)[2:], hex(self.face_type)[2:]))
            for fc in self.face_info:
                out_stream.write("\n{}".format(hex(fc[0])[2:]))
                cfn = len(fc)
                for cfi in range(1, cfn):
                    out_stream.write(" {}".format(hex(fc[cfi])[2:]))
            out_stream.write('))')
        else:
            out_stream.write("({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.face_type)[2:]))


@unique
class CellType(Enum):
    Dead = 0
    Fluid = 1
    Solid = 17


@unique
class CellElement(Enum):
    Mixed = 0
    Triangular = 1
    Tetrahedral = 2
    Quadrilateral = 3
    Hexahedral = 4
    Pyramid = 5
    Wedge = 6
    Polyhedral = 7


class XFCell(XFSection):
    def __init__(self, zone: int, first: int, last: int, ct: CellType, ce: CellElement, cell_info=None):
        """
        Cell Section in the ANSYS Fluent MSH File Format.
        :param zone: 区域号
        :param first: 起始序号(Starting from 1)
        :param last: 终止序号
        :param ct: Cell type.
        :param ce: Cell element.
        :param cell_info: 单元类型信息，一维数组，(last-first+1)个元素
        """

        super(XFCell, self).__init__(12)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.cell_type = ct.value
        self.elem_type = ce.value
        self.cell_info = None if cell_info is None else np.copy(cell_info)

    def build_content(self):
        if self.cell_info is not None:
            self.formatted_content = "({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:], hex(self.elem_type)[2:])
            self.formatted_content += "\n{}".format(hex(self.cell_info[0])[2:])
            for ci in range(1, len(self.cell_info)):
                self.formatted_content += " {}".format(hex(self.cell_info[ci])[2:])
            self.formatted_content += "))"
        else:
            if self.cell_type == 0:
                self.formatted_content = "({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:])  # declaration
            else:
                self.formatted_content = "({} ({} {} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:], hex(self.elem_type)[2:])

    @classmethod
    def declaration(cls, num):
        return cls(0, 1, num, CellType.Dead, CellElement.Mixed)

    def write(self, out_stream):
        if self.cell_info is not None:
            out_stream.write("({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:], hex(self.elem_type)[2:]))
            out_stream.write("\n{}".format(hex(self.cell_info[0])[2:]))
            cn = len(self.cell_info)
            for ci in range(1, cn):
                out_stream.write(" {}".format(hex(self.cell_info[ci])[2:]))
            out_stream.write("))")
        else:
            if self.cell_type == 0:
                out_stream.write("({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:]))  # declaration
            else:
                out_stream.write("({} ({} {} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:], hex(self.elem_type)[2:]))


class FluentMSH(object):
    def __init__(self):
        """
        ANSYS Fluent MSH File
        """

        self.xf_section = []

    def add(self, section):
        """
        向MSH文件添加信息段
        :param section: 信息段
        :type section: XFSection
        :return: None
        """

        self.xf_section.append(section)

    def clear(self):
        self.xf_section.clear()

    @property
    def size(self):
        return len(self.xf_section)

    def save(self, fn):
        """
        输出MSH文件
        :param fn: 输出文件名
        :type fn: str
        :return: None
        """

        f_out = open(fn, 'w')
        for i in range(self.size):
            self.xf_section[i].write(f_out)
            f_out.write('\n')
        f_out.close()

    @classmethod
    def from_nmf(cls, nmf):
        """
        Construct ANSYS Fluent MSH file from Neutral Map File.
        :param nmf: Neutral mapping description.
        :type nmf: NeutralMapFile
        :return: Structured grid in ANSYS Fluent MSH File Format.
        :rtype: FluentMSH
        """

        dim = nmf.dim
        cell_num = nmf.cell_num
        pnt_num = nmf.pnt_num
        face_num = nmf.face_num

        msh = cls()
        msh.add(XFHeader())
        msh.add(XFComment('Dimension:'))
        msh.add(XFDimension(dim))
        msh.add(XFComment("Cell Declaration:"))
        msh.add(XFCell.declaration(cell_num))
        msh.add(XFComment("Cell Info:"))
        zone_idx = 1
        msh.add(XFCell(zone_idx, 1, cell_num, CellType.Fluid, CellElement.Hexahedral))

        pnt_coordinate = np.empty((pnt_num, 3), float)
        cnt = 0
        for k in range(nmf.blk_num):
            u, v, w, _ = nmf.blk[k].shape
            cur_blk = nmf.blk[k]
            for ck in range(1, w-1):
                for cj in range(1, v-1):
                    for ci in range(1, u-1):
                        pnt_coordinate[cnt] = cur_blk[ci][cj][ck]
                        cnt += 1

        for k, blk in enumerate(blk_list):
            pn = shell_pnt_num[k]
            for pc in range(pn):
                if shell_pnt_idx[k][pc] != 0:
                    continue

                '''Record new pnt coordinate without duplication'''
                ci, cj, ck = shell_pnt_coord_from_idx(k, pc)
                dimensional_copy(pnt_desc[pnt_cnt], blk[ci][cj][ck], dim)
                pnt_cnt += 1

        return msh


class FluentMSHTestCase(unittest.TestCase):
    def test_2d(self):
        pass

    def test_3d(self):
        pass


if __name__ == '__main__':
    unittest.main()
