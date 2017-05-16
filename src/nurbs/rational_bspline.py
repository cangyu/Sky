import numpy as np
from src.nurbs.basis import FindSpan, BasisFuns
from scipy.misc import comb


def CurvePoint(n, p, U, Pw, u, C):
    """
    计算有理B样条曲线上的点
    :param n: 最后一个控制点下标
    :param p: 曲线的次数
    :param U: 节点序列
    :param Pw: 带权控制点序列
    :param u: 目标参数
    :param C: 曲线上u处的坐标
    :return: None
    """

    dim = len(Pw[0])
    N = np.zeros(p + 1)
    Cw = np.zeros(dim)

    span = FindSpan(n, p, u, U)
    BasisFuns(span, u, p, U, N)

    for j in range(0, p + 1):
        Cw += N[j] * Pw[span - p + j]

    for i in range(0, dim - 1):
        C[i] = Cw[i] / Cw[dim - 1]


def RatCurveDerive(Aders, wders, d, CK):
    """
    根据Cw(u)的导矢计算C(u)的导矢
    :param Aders: A(u)的0到d阶导矢，d+1个元素
    :param wders: w(u)的0到d阶导矢，d+1个元素
    :param d: 最大求导次数
    :param CK: 曲线上u处的点及其直到d阶的导矢，d+1个元素, CK[0]表示在u处的点
    :return: None
    """

    for k in range(0, d + 1):
        v = Aders[k]
        for i in range(1, k + 1):
            v -= comb(k, i) * wders[i] * CK[k - i]
        CK[k] = v / wders[0]


def SurfacePoint(n, p, U, m, q, V, Pw, u, v, S):
    """
    计算有理B样条曲面上的点
    :param n: u方向最后一个控制点下标
    :param p: u方向B样条基函数次数
    :param U: u方向节点矢量
    :param m: v方向最后一个控制点下标
    :param q: v方向B样条基函数次数
    :param V: v方向节点矢量
    :param Pw: 带权控制点，(n+1)*(m+1)个元素
    :param u: u方向参数
    :param v: v方向参数
    :param S: 曲面在(u,v)处的值
    :return: None
    """

    uspan = FindSpan(n, p, u, U)
    Nu = np.zeros(p + 1)
    BasisFuns(uspan, u, p, U, Nu)

    vspan = FindSpan(m, q, v, V)
    Nv = np.zeros(q + 1)
    BasisFuns(vspan, v, q, V, Nv)

    temp = np.zeros(q + 1, 4)
    for l in range(0, q + 1):
        for k in range(0, p + 1):
            temp[l] += Nu[k] * Pw[uspan - p + k][vspan - q + 1]

    Sw = np.zeros(4)
    for l in range(0, q + 1):
        Sw += Nv[l] * temp[l]

    S = np.zeros(3)
    for i in range(0, 3):
        S[i] = Sw[i] / Sw[3]


def RatSurfaceDerivs(Aders, wders, d, SKL):
    """
    计算曲面上S上(u,v)处的点及其直到d阶偏导矢
    :param Aders: 
    :param wders: 
    :param d: 
    :param SKL: 
    :return: 
    """

    for k in range(0, d + 1):
        for l in range(0, d + 1 - k):
            v = Aders[k][l]
            for j in range(1, l + 1):
                v -= comb(l, j) * wders[0][j] * SKL[k][l - j]
            for i in range(1, k + 1):
                v -= comb(k, i) * wders[i][0] * SKL[k - i][l]
                v2 = 0.0
                for j in range(1, l + 1):
                    v2 += comb(i, j) * wders[i][j] * SKL[k - i][l - j]
                v -= comb(k, i) * v2

            SKL[k][l] = v / wders[0][0]
