import numpy as np


def Horner1(a, n: int, u: float):
    """
    应用Horner规则计算幂基曲线上的点:
    C(u) = a0 + a1 * u + a2 * u^2 + ... + an * u^n
    :param a: 多项式系数列表，共有n+1个元素, 可多维
    :param n: 最后一个系数下标, 可为0
    :param u: 位置参数 
    :return: 在u处多项式的值
    """

    if len(a) != n + 1:
        raise ValueError("Invalid paramters!")

    ans = np.copy(a[n])
    i = n - 1
    while i >= 0:
        ans = ans * u + a[i]
        i -= 1

    return ans


def Bernstein(i: int, n: int, u: float):
    """
    以递推的方式计算Bernstein多项式的值:
    B(i,n,u) = C(i,n) * u^i * (1-u)^(n-i) = u * B(i-1, n-1, u) + (1-u) * B(i, n-1, u)
    :param i: 第i个Bernstein多项式
    :param n: Bernstein多项式次数, 可为0
    :param u: 无量纲位置参数
    :return: 第i个n次Bernstein多项式在u处的值
    """

    if u < 0.0 or u > 1.0:
        raise ValueError("Invalid parameter: \'u\' = " + str(u) + ", should be in range: [0, 1]")

    if i < 0 or i > n:
        return 0

    temp = np.zeros(n + 1, float)
    temp[n - i] = 1.0
    v = 1.0 - u
    k = 1
    while k <= n:
        j = n
        while j >= k:
            temp[j] = v * temp[j] + u * temp[j - 1]
            j -= 1
        k += 1

    return temp[n]


def AllBernstein(n: int, u: float, B):
    """
    计算所有n次Bernstein多项式的值
    :param n: Bernstein多项式次数, 可为0
    :param u: 无量纲位置参数
    :param B: 含有n+1个元素的一维数组:[B(0, n, u), B(1, n, u), B(2, n, u), ... , B(n, n, u)]
    :return: None
    """

    if u < 0.0 or u > 1.0:
        raise ValueError("Invalid parameter: \'u\' = " + str(u) + ", should be in range: [0, 1]")

    B[0] = 1.0
    v = 1.0 - u
    for j in range(1, n + 1):  # 从B(0, 0, u)=1 开始，迭代n层
        saved = 0
        for k in range(0, j):
            temp = B[k]
            B[k] = saved + v * temp
            saved = u * temp
        B[j] = saved


def PointOnBezierCurve(P, n, u, C):
    """
    计算Bezier曲线上的点
    :param P: Bezier曲线控制点序列, n+1个元素
    :param n: Bernstein多项式次数
    :param u: 无量纲位置参数
    :param C: 在u处点的坐标
    :return: None
    """

    if len(P) != n + 1:
        raise ValueError("Incompatible parameters!")

    B = np.zeros(n + 1, float)
    AllBernstein(n, u, B)

    C.fill(0)
    for k in range(0, n + 1):
        C += B[k] * P[k]


def deCasteljau1(P, n, u, C):
    """
    应用deCasteljau算法递推计算Bezier曲线上的点:
    C(P[0], P[1], ... , P[n]) = (1-u) * C(P[0], P[1], ... , P[n-1]) + u * C(P[1], P[2], ... , P[n])
    :param P: Bezier曲线控制点序列, n+1个元素
    :param n: Bernstein多项式次数
    :param u: 无量纲位置参数
    :param C: 在u处点的坐标(多维)
    :return: None
    """

    if len(P) != n + 1:
        raise ValueError("Incompatible parameters!")

    v = 1.0 - u
    Q = np.copy(P)
    for k in range(1, n + 1):
        for i in range(0, n - k + 1):
            Q[i] = v * Q[i] + u * Q[i + 1]

    C = np.copy(Q[0])


def Horner2(a, n: int, m: int, u: float, v: float):
    """
    计算幂基曲面上的点
    :param a: 幂基曲面多项式的常系数，(n+1)*(m+1)个元素, 可多维
    :param n: a的第1个dimension的最后一个元素下标，从0开始
    :param m: a的第2个dimension的最后一个元素下标，从0开始
    :param u: U方向的参数
    :param v: V方向的参数
    :return: (u, v)处的值
    """

    b = []
    for i in range(0, n + 1):
        b.append(Horner1(a[i], m, v))

    return Horner1(b, n, u)


def desCasteljau2(P, n: int, m: int, u: float, v: float):
    """
    利用deCasteljau算法递推计算Bezier曲面上的点
    :param P: 曲面的控制点，(n+1)*(m+1)个元素
    :param n: P的第1个dimension最后一个元素的下标，从0开始
    :param m: P的第2个dimension最后一个元素的下标，从0开始
    :param u: U方向的参数
    :param v: V方向的参数
    :return: (u, v)处的值
    """

    if u < 0 or u > 1.0:
        raise ValueError("Invalid parameter: \'u\' = " + str(u) + " should be within [0, 1]")
    if v < 0 or v > 1.0:
        raise ValueError("Invalid parameter: \'v\' = " + str(v) + " should be within [0, 1]")

    Bn = np.zeros(n + 1, float)
    Bm = np.zeros(m + 1, float)
    ans = np.copy(p[0][0])
    ans.fill(0.0)

    AllBernstein(n, u, Bn)
    AllBernstein(m, v, Bm)

    for i in range(0, n + 1):
        for j in range(0, m + 1):
            ans += Bm[j] * P[i][j]
        ans *= Bn[i]

    return ans
