import numpy as np
from scipy.optimize import newton


def single_exponential(N, A):
    """
    单指数分布
    :param N: 节点数量
    :param A: >0 : 节点向起始位置聚集
              <0 : 节点向终止位置聚集
    :return: [0,1]之间新的节点分布
    """

    if np.fabs(A) < 1e-12:
        raise AssertionError("\'A\' shouldn't be zero!")

    rho = np.linspace(0, 1.0, N)
    r = np.empty(N, float)

    ea1 = np.exp(A) - 1
    for i in range(N):
        r[i] = (np.exp(rho[i] * A) - 1) / ea1

    return r


def double_exponential(N, A1, A2, A3):
    """
    双指数分布，N个节点，曲线经过(A3, A1)
    :param N: 节点数量
    :param A1: 控制点纵坐标
    :param A2: >0 : 向起始点处聚集
               <0 : 向终止点处聚集
    :param A3: 控制点横坐标
    :return: [0,1]之间新的节点分布
    """

    p = (1 - A1) * A3 * (np.exp(A2) - 1) / ((1 - A3) * A1 * A2 * np.exp(A2))
    if p <= 0:
        raise AssertionError("Invalid parameters!")

    rho = np.linspace(0, 1.0, N)
    r = np.empty(N, float)

    p1z = np.log(p)  # 1阶导数零点
    A4 = newton(func=lambda x: np.exp(x) - 1.0 - p * x,
                x0=5 * p1z,
                fprime=lambda x: np.exp(x) - p,
                maxiter=20,
                fprime2=lambda x: np.exp(x))

    ea21 = np.exp(A2) - 1
    ea41 = np.exp(A4) - 1

    def f(x):
        if x <= A3:
            return A1 * (np.exp(A2 / A3 * x) - 1) / ea21
        else:
            return A1 + (1 - A1) * (np.exp(A4 / (1 - A3) * (x - A3)) - 1) / ea41

    for i in range(N):
        r[i] = f(rho[i])

    return r


def hyperbolic_tangent(N, B):
    """
    双曲正切分布
    :param N: 节点数量
    :param B: 控制参数
    :return: [0,1]之间新的节点分布
    """

    rho = np.linspace(0.0, 1.0, N)
    r = np.empty(N, float)

    tb = np.tanh(B)
    for i in range(N):
        r[i] = 1 + np.tanh(B * (rho[i] - 1)) / tb

    return r


def hyperbolic_sine(N, C):
    """
    双曲正弦分布
    :param N: 节点数量
    :param C: 控制参数
    :return: [0,1]之间新的节点分布
    """

    rho = np.linspace(0.0, 1.0, N)
    r = np.empty(N, float)

    sc = np.sinh(C)
    for i in range(N):
        r[i] = 1 + np.sinh(C * (1 - rho[i])) / sc

    return r
