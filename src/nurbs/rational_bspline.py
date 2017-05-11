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


