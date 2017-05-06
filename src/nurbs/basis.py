import numpy as np
import math

eps = 1e-6


def FindSpan(n: int, p: int, u: float, U):
    """
    确定参数u所在的节点区间的下标i, s.t. u belongs to [U[i], U[i+1])
    :param n: 最后一个控制点下标，从0开始
    :param p: NURBS曲线次数
    :param u: Target parameter
    :param U: 节点序列
    :return: u所在的节点区间下标i
    """

    # Corner case: u=U[m]，将其节点区间的下标设为n
    if u == U[n + 1]:
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
    计算在给定点u，所有p次非零B样条基函数的值(p+1个):$N_{i,p}(u)$
    :param i: 参数u所在的节点区间，s.t. u∈[U[i],U[i+1])
    :param u: Target parameter
    :param p: NURBS曲线次数
    :param U: 节点序列
    :param N: p+1个元素，因为只有 $N_{i-p,p}(u) ~ N_{i,p}(u)$ 不为0
    :return: None
    """

    left = np.zeros(p + 1, float)
    right = np.zeros(p + 1, float)

    N[0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - U[i + 1 - j]
        right[j] = U[i + j] - u
        saved = 0.0
        for r in range(0, j):
            tmp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * tmp
            saved = left[j - r] * tmp

        N[j] = saved


def DersBasisFuns(i: int, u: float, p: int, n: int, U, ders):
    """
    计算非零B样条基函数及其导数
    :param i: 
    :param u: 
    :param p: 
    :param n: 
    :param U: 
    :param ders: 
    :return: 
    """

    left = np.zeros(p + 1, float)
    right = np.zeros(p + 1, float)
    ndu = np.ones((p + 1, p + 1), float)

    for j in range(1, p + 1):
        left[j] = u - U[i + 1 - j]
        right[j] = U[i + j] - u
        saved = 0.0
        for r in range(0, j):
            # 下三角
            ndu[j][r] = right[r + 1] + left[j - r]
            temp = ndu[r][j - 1] / ndu[j][r]
            # 上三角
            ndu[r][j] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        ndu[j][j] = saved

    for j in range(0, p + 1):
        ders[0][j] = ndu[j][p]

    for r in range(0, p + 1):
        s1 = 0.0
        s2 = 1.0
        a[0][0] = 1.0
        for k in range(1, n + 1):
            d = 0.0
            rk = r - k
            pk = p - k

            if r >= k:
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk]
                d = a[s2][0] * ndu[rk][pk]

            if rk >= -1:
                j1 = 1
            else:
                j1 = -rk

            if r - 1 <= qk:
                j2 = k - 1
            else:
                j2 = p - r

            for j in range(j1, j2 + 1):
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j]
                d += a[s2][j] * ndu[rk + j][pk]

            if r <= pk:
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r]
                d += a[s2][k] * ndu[r][pk]

            ders[k][r] = d
            j = s1
            s1 = s2
            s2 = j

    r = p
    for k in range(1, n + 1):
        for j in range(0, p + 1):
            ders[k][j] *= r
        r *= (p - k)


def OneBasisFun(p, m, U, i, u, Nip):
    """
    计算基函数Nip
    :param p: 
    :param m: 
    :param U: 
    :param i: 
    :param u: 
    :param Nip: 
    :return: None
    """

    if (i == 0 and u == U[0]) or (i == m - p - 1 and u == U[m]):
        Nip = 1.0
        return

    if u < U[i] or u >= U[i + p + 1]:
        Nip = 0.0
        return

    for j in range(0, p + 1):
        if u >= U[i + j] and u < U[i + j + 1]:
            N[j] = 1.0
        else:
            N[j] = 0.0

    for k in range(1, p + 1):
        if math.fabs(N[0]) < eps:
            saved = 0.0
        else:
            saved = (u - U[i]) * N[0] / (U[i + k] - U[i])

        for j in range(0, p - k + 1):
            Uleft = U[i + j + 1]
            Uright = U[i + j + k + 1]
            if math.fabs(N[j + 1]) < eps:
                N[j] = saved
                saved = 0.0
            else:
                temp = N[j + 1] / (Uright - Uleft)
                N[j] = saved + (Uright - u) * temp
                saved = (u - Uleft) * temp
    Nip = N[0]


def DersOneBasisFun(p, m, U, i, u, n, ders):
    """
    计算基函数Nip的各阶导数
    :param p: 
    :param m: 
    :param U: 
    :param i: 
    :param u: 
    :param n: 
    :param ders: 
    :return: 
    """

    if u < U[i] or u >= U[i + p + 1]:
        for k in range(0, n + 1):
            ders[k] = 0.0
        return

    for j in range(0, p + 1):
        if u >= U[i + j] and u < U[i + j + 1]:
            N[j][0] = 1.0
        else:
            N[j][0] = 0.0

    for k in range(1, p + 1):
        if math.fabs(N[0][k - 1]) < eps:
            saved = 0.0
        else:
            saved = (u - U[i]) * N[0][k - 1] / (U[i + k] - U[i])

        for j in range(0, p - k + 1):
            Uleft = U[i + j + 1]
            Uright = U[i + j + k + 1]
            if math.fabs(N[j + 1][k - 1]) < eps:
                N[j][k] = saved
                saved = 0.0
            else:
                temp = N[j + 1][k - 1] / (Uright - Uleft)
                N[j][k] = saved + (Uright - u) * temp
                saved = (u - Uleft) * temp

    ders[0] = N[0][p]
    for k in range(1, n + 1):
        for j in range(0, k + 1):
            ND[j] = N[j][p - k]
        for jj in range(1, k + 1):
            if math.fabs(ND[0]) < eps:
                saved = 0.0
            else:
                saved = ND[0] / (U[i + p - k + jj] - U[i])

            for j in range(0, k - jj + 1):
                Uleft = U[i + j + 1]
                Uright = U[i + j + p + jj + 1]
                if math.fabs(ND[j + 1]) < eps:
                    ND[j] = (p - k + jj) * saved
                    saved = 0.0
                else:
                    temp = ND[j + 1] / (Uright - Uleft)
                    ND[j] = (p - k + jj) * (saved - temp)
                    saved = temp
        ders[k] = ND[0]
