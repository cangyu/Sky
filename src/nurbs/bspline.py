import numpy as np
from src.nurbs.basis import FindSpan, BasisFuns, AllBasisFuns, DersBasisFuns


def CurvePoint(n: int, p: int, U, P, u: float, C):
    """
    计算B样条曲线上的点
    :param n: 最后一个控制点下标
    :param p: B样条曲线的次数 
    :param U: 节点序列
    :param P: 控制点序列,(n+1)个元素
    :param u: 目标参数
    :param C: 曲线在u处的值
    :return: None
    """

    N = np.zeros((p + 1, len(C)))
    span = FindSpan(n, p, u, U)
    BasisFuns(span, u, p, U, N)
    C.fill(0)
    for i in range(0, p + 1):
        C += N[i] * P[span - p + i]


def CurveDerivsAlg1(n: int, p: int, U, P, u: float, d: int, CK):
    """
    计算曲线上的导矢
    :param n: 最后一个控制点下标
    :param p: B样条曲线的次数
    :param U: 节点序列
    :param P: 控制点序列,(n+1)个元素
    :param u: 目标参数
    :param d: 计算直到d阶导矢
    :param CK: d+1个元素，CK[k]返回k阶导矢
    :return: None
    """

    du = min(d, p)
    nders = np.zeros((du + 1, p + 1))

    for k in range(p + 1, d + 1):
        CK[k] = 0.0

    span = FindSpan(n, p, u, U)
    DersBasisFuns(span, u, p, du, U, nders)

    for k in range(0, du + 1):
        CK[k] = 0.0
        for j in range(0, p + 1):
            CK[k] += nders[k][j] * P[span - p + j]


def CurveDerivCpts(n, p, U, P, d, r1, r2, PK):
    """
    计算曲线直到d阶的所有导曲线的控制点
    :param n: 最后一个控制点下标
    :param p: B样条曲线的次数
    :param U: 节点序列
    :param P: 控制点序列,(n+1)个元素
    :param d: 最大求导阶数(d<=p)
    :param r1: 待求控制点的起始下标
    :param r2: 待求控制点的终止下标
    :param PK: (d+1)*(r2-r1+1)个元素，PK[k][i]返回k阶导曲线的第i个控制点
    :return: None
    """

    r = r2 - r1
    for i in range(0, r + 1):
        PK[0][i] = P[r1 + i]

    for k in range(1, d + 1):
        tmp = p - k + 1
        for i in range(0, r - k + 1):
            PK[k][i] = tmp * (PK[k - 1][i + 1] - PK[k - 1][i]) / (U[r1 + i + p + 1] - U[r1 + i + k])


def CurveDerivsAlg2(n, p, U, P, u, d, CK):
    """
    计算B样条曲线上的点直到d阶导矢，和CurveDerivsAlg1目的相同，但用的公式不同
    :param n: 最后一个控制点下标
    :param p: B样条曲线的次数
    :param U: 节点序列
    :param P: 控制点序列,n+1个元素
    :param u: 目标参数
    :param d: 最大求导阶数(d<=p)
    :param CK: d+1个元素，B样条曲线上u处的点及其直到d阶的导矢
    :return: None
    """

    dim = len(p[0])
    du = min(d, p)
    N = np.zeros((p + 1, p + 1))
    PK = np.zeros((du + 1, p + 1, dim))

    for k in range(p + 1, d + 1):
        CK[k].fill(0.0)

    span = FindSpan(n, p, u, U)
    AllBasisFuns(span, u, p, U, N)
    CurveDerivCpts(n, p, U, P, du, span - p, span, PK)

    for k in range(0, du + 1):
        CK[k].fill(0.0)
        for j in range(0, p - k + 1):
            CK[k] += N[j][p - k] * PK[k][j]


def SurfacePoint(n, p, U, m, q, V, P, u, v, S):
    """
    计算B样条曲面上的点
    :param n: u方向最后一个控制点下标
    :param p: u方向B样条基函数次数
    :param U: u方向节点序列
    :param m: v方向最后一个控制点下标
    :param q: v方向B样条基函数次数
    :param V: v方向节点序列
    :param P: 控制点序列，(n+1)*(m+1)个元素
    :param u: u方向参数
    :param v: v方向参数
    :param S: (u,v)处坐标
    :return: None
    """

    Nu = np.zeros(p + 1)
    Nv = np.zeros(q + 1)

    uspan = FindSpan(n, p, u, U)
    BasisFuns(uspan, u, p, U, Nu)
    vspan = FindSpan(m, q, v, V)
    BasisFuns(vspan, v, q, V, Nv)
    uind = uspan - p

    temp = np.zeros((q + 1, len(P[0][0])))
    for l in range(0, q + 1):
        vind = vspan - q + 1
        for k in range(0, p + 1):
            temp[l] += Nu[k] * P[uind + k][vind]

    S.fill(0.0)
    for l in range(0, q + 1):
        S += Nv[l] * temp[l]


def SurfaveDerivsAlg1(n, p, U, m, q, V, P, u, v, d, SKL):
    """
    计算B样条曲面上的点及其所有直到d阶的偏导矢
    :param n: u方向最后一个控制点下标
    :param p: u方向B样条基函数次数
    :param U: u方向节点序列
    :param m: v方向最后一个控制点下标
    :param q: v方向B样条基函数次数
    :param V: v方向节点序列
    :param P: 控制点序列，(n+1)*(m+1)个元素
    :param u: u方向参数
    :param v: v方向参数
    :param d: 最大求导次数
    :param SKL: 所有直到d阶的偏导矢，(d+1)*(d+1)个元素
    :return: None
    """

    du = min(d, p)
    dv = min(d, q)

    SKL.fill(0.0)

    Nu = np.zeros((du + 1, p + 1))
    Nv = np.zeros((dv + 1, q + 1))

    uspan = FindSpan(n, p, u, U)
    DersBasisFuns(uspan, u, p, du, U, Nu)
    vspan = FindSpan(m, q, v, V)
    DersBasisFuns(vspan, v, q, dv, V, Nv)

    dim = len(P[0][0])
    temp = np.zeros((q + 1, dim))
    for k in range(0, du + 1):
        for s in range(0, q + 1):
            for r in range(0, p + 1):
                temp[s] += Nu[k][r] * P[uspan - p + r][vspan - q + s]

        dd = min(d - k, dv)
        for l in range(0, dd + 1):
            for s in range(0, q + 1):
                SKL[k][l] += Nv[l][s] * temp[s]


def SurfaceDerivCpts(n, p, U, m, q, V, P, d, r1, r2, s1, s2, PKL):
    """
    计算导曲面的控制点
    :param n: 
    :param p: 
    :param U: 
    :param m: 
    :param q: 
    :param V: 
    :param P: 
    :param d: 
    :param r1: 
    :param r2: 
    :param s1: 
    :param s2: 
    :param PKL: 
    :return: 
    """

    du = min(d, p)
    dv = min(d, q)

    r = r2 - r1
    s = s2 - s1

    for j in range(s1, s2 + 1):
        CurveDerivCpts(n, p, U, & P[][j], du, r1, r2, temp)
        for k in range(0, du + 1):
            for i in range(0, r - k + 1):
                PKL[k][0][i][j - s1] = temp[k][l]

    for k in range(0, du):
        for i in range(0, r - k + 1):
            dd = min(d - k, dv)
            CurveDerivCpts(m, q, & V[s1], PKL[k][0][i], dd, 0, s, temp)
            for l in range(1, dd + 1):
                for j in range(0, s - l + 1):
                    PKL[k][k][i][j] = temp[l][j]


def SurfaceDerivsAlg2(n, p, U, m, q, V, P, u, v, d, SKL):
    """
    计算B样条曲面上的点及所有直到d阶的偏导矢
    :param n: 
    :param p: 
    :param U: 
    :param m: 
    :param q: 
    :param V: 
    :param P: 
    :param u: 
    :param v: 
    :param d: 
    :param SKL: 
    :return: 
    """

    du = min(d, p)
    for k in range(p + 1, d + 1):
        for l in range(0, d - k + 1):
            SKL[k][l] = 0.0

    dv = min(d, q)
    for k in range(q + 1, d + 1):
        for l in range(0, d - l + 1):
            SKL[k][l] = 0.0

    uspan = FindSpan(n, p, u, U)
    AllBernstein(uspan, u, p, U, Nu)
    vspan = FindSpan(m, q, v, V)
    AllBernstein(vspan, v, q, v, Nv)
    SurfaceDerivCpts(n, p, U, m, q, V, P, d, uspan - p, uspan, vspan - q, vspan, PKL)
    for k in range(0, du + 1):
        dd = min(d - k, dv)
        for l in range(0, dd + 1):
            SKL[k][l] = 0.0
            for i in range(0, q - l + 1):
                tmp = 0.0
                for j in range(0, p - k + 1):
                    tmp += Nu[j][p - k] * PKL[k][l][j][i]
                SKL[k][l] += Nv[i][q - 1] * tmp
