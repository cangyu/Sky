import numpy as np
import math

EPS = 1e-9


def equal(a: float, b: float):
    """
    判断两个浮点数是否相等
    """

    if math.fabs(a - b) < EPS:
        return True
    else:
        return False


def pnt_dist(lhs, rhs):
    """
    计算两点之间距离
    """

    dim = min(len(lhs), len(rhs))
    ans = 0.0
    for i in range(0, dim):
        ans += math.pow(lhs[i] - rhs[i], 2)

    ans = math.sqrt(ans)
    return ans


def to_homogeneous(pnt, w=1.0):
    """
    转换为齐次坐标
    :param pnt: 原始坐标
    :param w: 权系数
    :return: 带权齐次坐标
    """

    pw = np.zeros(len(pnt) + 1)
    pw[-1] = w
    if equal(w, 0.0):
        for i in range(0, len(pnt)):
            pw[i] = pnt[i]
    else:
        for i in range(0, len(pnt)):
            pw[i] = pnt[i] * w

    return pw


def to_cartesian(pnt):
    """
    转换为笛卡尔坐标
    :param pnt: 原始齐次坐标
    :return: 投影后得到的笛卡尔坐标
    """

    p = np.zeros(len(pnt) - 1)
    if equal(pnt[-1], 0.0):
        for i in range(0, len(p)):
            p[i] = pnt[i]
    else:
        for i in range(0, len(p)):
            p[i] = pnt[i] / pnt[-1]

    return p


def find_span(n: int, p: int, u: float, U):
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
        mid = int((high + low) / 2)

    return mid


def all_basis_val(u, p, U):
    """
    计算在给定点u，所有p次非零B样条基函数的值:N(i, p, u)
    :param u: 目标参数
    :param p: B样条基函数的次数
    :param U: 节点序列
    :return: n+1个元素，因为只有N(i-p,p,u) ~ N(i,p,u)不为0
    """

    m = len(U) - 1
    n = m - p - 1
    i = find_span(n, p, u, U)  # 参数u所在的节点区间，s.t. u∈[U[i],U[i+1])

    N = np.zeros(p + 1, float)
    left = np.zeros(p + 1, float)
    right = np.zeros(p + 1, float)

    N[0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - U[i + 1 - j]
        right[j] = U[i + j] - u
        saved = 0.0
        for r in range(0, j):
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        N[j] = saved

    m = len(U) - 1
    n = m - p - 1
    ans = np.zeros(n + 1, float)
    for k in range(i - p, i + 1):
        ans[k] = N[k - (i - p)]

    return ans
