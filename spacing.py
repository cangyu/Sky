import unittest
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy
from scipy.optimize import newton


def linear_expand(seq, begin, end):
    """
    将序列线性缩放到给定的起始值
    :param seq: 待缩放序列
    :param begin: 缩放后的起点
    :type begin: float
    :param end: 缩放后的终点
    :type end:float
    :return: 缩放后的序列
    """

    r_min = min(seq)
    r_max = max(seq)
    dlt = end - begin
    r_dlt = r_max - r_min
    return np.full_like(seq, begin) if r_dlt == 0 else np.copy(list(map(lambda u: (u - r_min) / r_dlt * dlt + begin, seq)))


def uniform(n):
    """
    均匀分布
    :param n: 节点数量
    :type n: int
    :return: [0,1]之间均匀分布序列
    """

    return np.linspace(0, 1, n)


def chebshev_dist(start, end, n):
    """
    生成切比雪夫点
    :param start: 起始值
    :type start: float
    :param end: 终止值
    :type end: float
    :param n: 采样点数量
    :type n: int
    :return: n个点
    """

    ang = np.linspace(np.pi, 0, n)
    pr = np.cos(ang)
    for i in range(0, n):
        pr[i] = start + (end - start) / 2 * (pr[i] + 1)

    return pr


def chebshev_dist_multi(seg, num):
    """
    多段Chebshev分布
    :param seg: 分段点
    :param num: 每个分段内点的数量(包括首尾)
    :return: 多段Chebshev分布数组
    """

    if len(seg) != len(num) + 1:
        raise AssertionError("Unmatched settings.")

    ret = chebshev_dist(seg[0], seg[1], num[0])
    for k in range(1, len(num)):
        csg = chebshev_dist(seg[k], seg[k + 1], num[k])
        ret = np.concatenate((ret, csg[1:]))

    return ret


def single_exponential(n, pa):
    """
    单指数分布
    :param n: 节点数量
    :type n: int
    :param pa: >0 节点向起始位置聚集, <0 节点向终止位置聚集
    :type pa: float
    :return: [0,1]之间新的单指数分布
    """

    if math.isclose(pa, 0):
        raise AssertionError("\'pa\' shouldn't be zero!")

    rho = np.linspace(0, 1.0, n)
    r = np.empty(n, float)

    ea1 = np.exp(pa) - 1
    for i in range(n):
        r[i] = (np.exp(rho[i] * pa) - 1) / ea1

    return r


def double_exponential(n, pa1, pa2, pa3):
    """
    双指数分布，n个节点，曲线经过(pa3, pa1)
    :param n: 节点数量
    :type n: int
    :param pa1: 控制点纵坐标
    :type pa1: float
    :param pa2: >0 向起始点处聚集, <0 向终止点处聚集
    :type pa2: float
    :param pa3: 控制点横坐标
    :type pa3: float
    :return: [0,1]之间的双指数分布
    """

    p = (1 - pa1) * pa3 * (np.exp(pa2) - 1) / ((1 - pa3) * pa1 * pa2 * np.exp(pa2))
    if p <= 0:
        raise AssertionError("Invalid parameters!")

    rho = np.linspace(0, 1.0, n)
    r = np.empty(n, float)

    p1z = np.log(p)  # 1阶导数零点
    pa4 = newton(func=lambda x: np.exp(x) - 1.0 - p * x,
                 x0=5 * p1z,
                 fprime=lambda x: np.exp(x) - p,
                 maxiter=50,
                 fprime2=lambda x: np.exp(x))

    ea21 = np.exp(pa2) - 1
    ea41 = np.exp(pa4) - 1

    def f(x):
        if x <= pa3:
            return pa1 * (np.exp(pa2 / pa3 * x) - 1) / ea21
        else:
            return pa1 + (1 - pa1) * (np.exp(pa4 / (1 - pa3) * (x - pa3)) - 1) / ea41

    for i in range(n):
        r[i] = f(rho[i])

    return r


def hyperbolic_tangent(n, pb):
    """
    双曲正切分布
    :param n: 节点数量
    :type n: int
    :param pb: 控制参数
    :type pb: float
    :return: [0,1]之间的双曲正切分布
    """

    rho = np.linspace(0.0, 1.0, n)
    r = np.empty(n, float)

    tb = np.tanh(pb)
    for i in range(n):
        r[i] = 1 + np.tanh(pb * (rho[i] - 1)) / tb

    return r


def hyperbolic_sine(n, pc):
    """
    双曲正弦分布
    :param n: 节点数量
    :type n: int
    :param pc: 控制参数
    :type pc: float
    :return: [0,1]之间的双曲正弦分布
    """

    rho = np.linspace(0.0, 1.0, n)
    r = np.empty(n, float)

    sc = np.sinh(pc)
    for i in range(n):
        r[i] = 1 + np.sinh(pc * (rho[i] - 1)) / sc

    return r


class SpacingTester(unittest.TestCase):
    def test_single_exponential(self):
        # n, pa
        data = [(101, 6),
                (101, 5),
                (101, 4)]
        ans = []

        for i in range(len(data)):
            r = single_exponential(data[i][0], data[i][1])
            plt.figure()
            plt.plot(np.linspace(0.0, 1.0, data[i][0]), r)
            plt.xlabel('u')
            plt.ylabel('f(u)')
            plt.title('Single Exponential\nA={}'.format(data[i][1]))
            plt.show()
            ans.append(plt)

        self.assertTrue(len(data) == len(ans))

    def test_double_exponential(self):
        # n, pa1, pa2, pa3
        data = [(401, 0.5, -1.5, 0.5),  # A2取负值使得中间较密
                (401, 0.5, 1.5, 0.5)]  # A2取正值使得中间稀疏，两边密集
        ans = []

        for i in range(len(data)):
            r = double_exponential(data[i][0], data[i][1], data[i][2], data[i][3])
            plt.figure()
            plt.plot(np.linspace(0.0, 1.0, data[i][0]), r)
            plt.xlabel('u')
            plt.ylabel('f(u)')
            plt.title('Double Exponential\nA1={}, A2={}, A3={}'.format(data[i][1], data[i][2], data[i][3]))
            plt.show()
            ans.append(plt)

        self.assertTrue(len(data) == len(ans))

    def test_hyperbolic_tangent(self):
        # n, pb
        data = [(201, 3),
                (201, 2.5)]
        ans = []

        for k, dt in enumerate(data):
            r = hyperbolic_tangent(dt[0], dt[1])
            plt.figure()
            plt.plot(np.linspace(0.0, 1.0, dt[0]), r)
            plt.xlabel('u')
            plt.ylabel('f(u)')
            plt.title('Hyperbolic Tangent\nB={}'.format(dt[1]))
            plt.show()
            ans.append(plt)

        self.assertTrue(len(data) == len(ans))

    def test_hyperbolic_sine(self):
        # n, pc
        data = [(201, 2),
                (201, 2.5),
                (201, 3)]
        ans = []

        for k, dt in enumerate(data):
            r = hyperbolic_sine(dt[0], dt[1])
            plt.figure()
            plt.plot(np.linspace(0.0, 1.0, dt[0]), r)
            plt.xlabel('u')
            plt.ylabel('f(u)')
            plt.title('Hyperbolic Sine\nC={}'.format(dt[1]))
            plt.show()
            ans.append(plt)

        self.assertTrue(len(data) == len(ans))

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


if __name__ == '__main__':
    unittest.main()
