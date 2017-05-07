import numpy as np
import math
from .basic import AllBernstein
from .basis import FindSpan, BasisFuns, DersBasisFuns


def CurvePoint(n, p, U, P, u, C):
    """
    计算B样条曲线上的点
    :param n: 
    :param o: 
    :param U: 
    :param P: 
    :param u: 
    :param C: 
    :return: 
    """

    N = np.zeros((p + 1, len(C)))
    span = FindSpan(n, p, u, U)
    BasisFuns(span, u, p, U, N)
    C.fill(0)
    for i in range(0, p + 1):
        C += N[i] * P[span - p + i]


def CurveDerivsAlg1(n, p, U, P, u, d, CK):
    """
    计算曲线上的导矢
    :param n: 
    :param p: 
    :param U: 
    :param P: 
    :param u: 
    :param d: 
    :param CK: 
    :return: 
    """

    du = min(d, p)
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
    计算导曲线的控制点
    :param n: 
    :param p: 
    :param U: 
    :param P: 
    :param d: 
    :param r1: 
    :param r2: 
    :param PK: 
    :return: 
    """

    r = r2 - r1
    for i in range(0, r + 1):
        PK[0][i] = P[r1 + i]

    for k in range(1, d + 1):
        tmp = p - k + 1
        for i in range(0, r - k + 1):
            PK[k][i] = tmp * (PK[k - 1][i + 1] - PK[k - 1][i]) / (U[r1 + i + p + 1] - U[r1 + i + k])


def CurveDerivsAlg2(u, p, U, P, d, u, d, CK):
    """
    计算B样条曲线上的点直到d阶导矢
    :param u: 
    :param p: 
    :param U: 
    :param P: 
    :param d: 
    :param CK: 
    :return: 
    """

    du = min(d, p)
    for k in range(p + 1, d + 1):
        CK[k] = 0.0

    span = FindSpan(n, p, u, U)
    AllBernstein(span, u, p, U, N)
    CurveDerivCpts(n, p, U, P, du, span - p, span, PK)
    for k in range(0, du + 1):
        CK[k] = 0.0
        for j in range(0, p - k + 1):
            CK[k] += N[j][p - k] * PK[k][j]


def SurfacePoint(n, p, U, m, q, V, P, u, v, S):
    """
    计算B样条曲面上的点
    :param n: 
    :param p: 
    :param U: 
    :param m: 
    :param q: 
    :param V: 
    :param P: 
    :param u: 
    :param v: 
    :param S: 
    :return: 
    """

    uspan = FindSpan(n, p, u, U)
    BasisFuns(uspan, u, p, U, Nu)
    vspan = FindSpan(m, q, v, V)
    BasisFuns(vspan, v, q, V, Nv)

    for l in range(0, q + 1):
        temp[l] = 0.0
        vind = vspan - q + 1
        for k in range(0, p + 1):
            temp[l] += Nu[k] * P[uind + k][vind]

    S = 0.0
    for l in range(0, q + 1):
        S += Nv[l] * temp[l]


def SurfaveDerivsAlg1(n, p, U, m, q, V, P, u, v, d, SKL):
    """
    计算B样条曲面上的点及其所有直到d阶的偏导矢
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
    for l in range(q + 1, d + 1):
        for k in range(0, d - l + 1):
            SKL[k][l] = 0.0

    uspan = FindSpan(n, p, u, U)
    DersBasisFuns(uspan, u, p, du, U, Nu)
    vspan = FindSpan(m, q, v, V)
    DersBasisFuns(vspan, v, q, dv, V, Nv)

    for k in range(0, du + 1):
        for s in range(0, q + 1):
            temp[s] = 0.0
            for r in range(0, p + 1):
                temp[s] += Nu[k][r] * P[uspan - p + r][vspan - q + s]

        dd = min(d - k, dv)
        for l in range(0, dd + 1):
            SKL[k][l] = 0.0
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
