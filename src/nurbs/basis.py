import numpy as np
import math

eps = 1e-6


def equal(a: float, b: float):
    return True if math.fabs(a - b) < eps else False


def FindSpan(n: int, p: int, u: float, U):
    """
    确定参数u所在的节点区间的下标i, s.t. u belongs to [U[i], U[i+1])
    :param n: 最后一个控制点下标，从0开始
    :param p: B样条基函数的次数
    :param u: Target parameter
    :param U: 节点序列
    :return: u所在的节点区间下标i
    """

    if u == U[n + 1]:  # Corner case: u=U[m]，将其节点区间的下标设为n
        return n

    low = p
    high = n + 1
    mid = int((high + low) / 2)
    while u < U[mid] or u >= U[mid + 1]:
        if u < U[mid]:
            high = mid
        else:
            low = mid
        mid = (low + high) / 2

    return mid


def BasisFuns(i: int, u: float, p: int, U, N):
    """
    计算在给定点u，所有p次非零B样条基函数的值:N(i, p, u)
    :param i: 参数u所在的节点区间，s.t. u∈[U[i],U[i+1])
    :param u: 目标参数
    :param p: B样条基函数的次数
    :param U: 节点序列
    :param N: p+1个元素，因为只有N(i-p,p,u) ~ N(i,p,u)不为0
    :return: None
    """

    left = np.zeros(p + 1, float)
    right = np.zeros(p + 1, float)

    N[0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - U[i + 1 - j]  # left[k], right[k]仅在计算第k阶基函数时用到，所以可以逐次计算
        right[j] = U[i + j] - u
        saved = 0.0
        for r in range(0, j):
            tmp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * tmp
            saved = left[j - r] * tmp

        N[j] = saved


def AllBasisFuns(i, u, p, U, N):
    """
    计算所有0到p次非0基函数的值
    :param i: 参数u所在的节点区间，s.t. u∈[U[i],U[i+1])
    :param u: 目标参数
    :param p: 最大B样条基函数次数
    :param U: 节点序列
    :param N: 所有0到p次非0基函数的值,(p+1)*(p+1)个元素，N[jj][ii]存储ii次基函数第jj个非0值，ii<=p, jj<=ii
    :return: None
    """

    for pp in range(0, p + 1):
        '''计算次数为pp时所有非零B样条基函数的值'''
        data = np.zeros(pp + 1)
        BasisFuns(i, u, pp, U, data)

        '''Copy result'''
        for k in range(0, pp + 1):
            N[k][pp] = data[k]


def DersBasisFuns(i: int, u: float, p: int, n: int, U, ders):
    """
    计算非零B样条基函数及其导数
    :param i: 参数u所在的节点区间，s.t. u∈[U[i],U[i+1])
    :param u: 目标参数
    :param p: B样条基函数的次数
    :param n: 输出非零基函数及其直到n阶导数(n<=p)
    :param U: 节点序列
    :param ders: (n+1)*(p+1)个元素，ders[k][j] 存储N(i-p+j, p)的k阶导数
    :return: None
    """

    left = np.zeros(p + 1)
    right = np.zeros(p + 1)
    ndu = np.zeros((p + 1, p + 1))  # 存储基函数与节点之差
    a = np.zeros((2, p + 1))  # 以不断更新的方式存储最新的两行a[k][j], a[k-1][j]

    '''将函数值与节点差存放到数组中'''
    ndu[0][0] = 1.0
    for j in range(1, p + 1):  # 逐个矩形层推进
        left[j] = u - U[i + 1 - j]
        right[j] = U[i + j] - u
        saved = 0.0
        for r in range(0, j):
            '''下三角'''
            ndu[j][r] = right[r + 1] + left[j - r]
            temp = ndu[r][j - 1] / ndu[j][r]
            '''上三角'''
            ndu[r][j] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        ndu[j][j] = saved

    '''载入基函数的值'''
    for j in range(0, p + 1):
        ders[0][j] = ndu[j][p]

    '''计算导数'''
    for r in range(0, p + 1):  # 对函数下标进行循环
        s1 = 0
        s2 = 1
        a[0][0] = 1.0
        for k in range(1, n + 1):  # 循环计算k阶导数
            d = 0.0
            rk = r - k
            pk = p - k

            if r >= k:
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk]
                d = a[s2][0] * ndu[rk][pk]

            j1 = 1 if rk >= -1 else -rk
            j2 = k - 1 if r - 1 <= pk else p - r

            for j in range(j1, j2 + 1):
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j]
                d += a[s2][j] * ndu[rk + j][pk]

            if r <= pk:
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r]
                d += a[s2][k] * ndu[r][pk]

            ders[k][r] = d

            '''转换行'''
            j = s1
            s1 = s2
            s2 = j

    '''乘以正确的因子'''
    r = p
    for k in range(1, n + 1):
        for j in range(0, p + 1):
            ders[k][j] *= r
        r *= (p - k)


def OneBasisFun(p: int, m: int, U, i: int, u: float):
    """
    计算基函数N(i, p)在u处的值
    :param p: B样条基函数的次数
    :param m: 最后一个节点的下标
    :param U: 节点序列，m+1个元素
    :param i: 参数u所在的节点区间，s.t. u∈[U[i],U[i+1])
    :param u: 目标参数
    :return: 基函数N(i, p)在u处的值
    """

    N = np.zeros(p + 1)

    if (i == 0 and equal(u, U[0])) or (i == m - p - 1 and equal(u, U[m])):  # 特殊情况
        return 1.0

    if u < U[i] or u >= U[i + p + 1]:  # 局部支撑
        return 0.0

    '''初始化0次基函数'''
    for j in range(0, p + 1):
        if u >= U[i + j] and u < U[i + j + 1]:
            N[j] = 1.0

    '''计算三角形表'''
    for k in range(1, p + 1):
        saved = 0.0 if equal(N[0], 0.0) else (u - U[i]) * N[0] / (U[i + k] - U[i])

        for j in range(0, p - k + 1):
            Uleft = U[i + j + 1]
            Uright = U[i + j + k + 1]
            if equal(N[j + 1], 0.0):
                N[j] = saved
                saved = 0.0
            else:
                temp = N[j + 1] / (Uright - Uleft)
                N[j] = saved + (Uright - u) * temp
                saved = (u - Uleft) * temp

    return N[0]


def DersOneBasisFun(p, m, U, i, u, n, ders):
    """
    计算基函数Nip的各阶导数， 先计算出0次基函数，再由低次基函数导数值递推高次基函数导数值
    :param p: B样条基函数的次数
    :param m: 最后一个节点的下标
    :param U: 节点序列，m+1个元素
    :param i: 参数u所在的节点区间，s.t. u∈[U[i],U[i+1])
    :param u: 目标参数
    :param n: 输出非零基函数及其直到n阶导数(n<=p)
    :param ders: n+1个元素，ders[k]返回k阶导数
    :return: None
    """

    N = np.zeros((p + 1, p + 1))
    ND = np.zeros(n + 1)

    if u < U[i] or u >= U[i + p + 1]:  # 局部支撑性
        for k in range(0, n + 1):
            ders[k] = 0.0
        return

    '''初始化0次基函数'''
    for j in range(0, p + 1):
        if u >= U[i + j] and u < U[i + j + 1]:
            N[j][0] = 1.0

    '''计算完整的三角形表'''
    for k in range(1, p + 1):
        saved = 0.0 if equal(N[0][k - 1], 0.0) else (u - U[i]) * N[0][k - 1] / (U[i + k] - U[i])
        for j in range(0, p - k + 1):
            Uleft = U[i + j + 1]
            Uright = U[i + j + k + 1]
            if equal(N[j + 1][k - 1], 0.0):
                N[j][k] = saved
                saved = 0.0
            else:
                temp = N[j + 1][k - 1] / (Uright - Uleft)
                N[j][k] = saved + (Uright - u) * temp
                saved = (u - Uleft) * temp

    '''函数值'''
    ders[0] = N[0][p]

    '''计算导数'''
    for k in range(1, n + 1):
        '''载入正确的列'''
        for j in range(0, k + 1):
            ND[j] = N[j][p - k]

        '''计算宽度为k的表'''
        for jj in range(1, k + 1):
            saved = 0.0 if equal(ND[0], 0.0) else ND[0] / (U[i + p - k + jj] - U[i])
            for j in range(0, k - jj + 1):
                Uleft = U[i + j + 1]
                Uright = U[i + j + p + jj + 1]
                if equal(ND[j + 1], 0.0):
                    ND[j] = (p - k + jj) * saved
                    saved = 0.0
                else:
                    temp = ND[j + 1] / (Uright - Uleft)
                    ND[j] = (p - k + jj) * (saved - temp)
                    saved = temp

        ders[k] = ND[0]  # k阶导数
