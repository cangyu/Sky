import numpy as np
import math

EPS = 1e-9


def equal(a: float, b: float):
    """
    判断两个浮点数是否相等
    """

    return math.fabs(a - b) < EPS


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


def array_smart_copy(src, dst):
    for i in range(min(len(src), len(dst))):
        dst[i] = src[i]


def line_intersection(p1, u, p2, v, with_ratio=False):
    """
    Calculate the intersection of 2 straight lines.
    :param p1: Point on the first line.
    :param u: Direction vector of the first line.
    :param p2: Point on the second line.
    :param v: Direction vector of the second line.
    :param with_ratio: To indicate if additional parameters are returned or not.
    :return:
        If with_ratio=True, then return (alpha1, alpha2, P), such that P = p1 + alpha1 * u = p2 + alpha2 * v,
        otherwise, only the intersection point P is returned.
    """

    dim = len(p1)
    if len(p1) != len(u) != len(p2) != len(v):
        raise ValueError("Inconsistent dimension!")
    if dim not in (2, 3):
        raise ValueError("Invalid dimension!")

    P1 = np.copy(p1)
    U = np.copy(u)
    P2 = np.copy(p2)
    V = np.copy(v)
    dP = P2 - P1

    if not U.any():
        raise AssertionError("Invalid U direction vector!")
    if not V.any():
        raise AssertionError("Invalid V direction vector!")
    if dim == 3 and np.inner(dP, np.cross(U, V + dP)).any():
        raise AssertionError("Two lines are non-coplanar!")
    if not np.cross(U, V).any():
        raise AssertionError("No intersection!")

    U90 = np.zeros(dim)
    V90 = np.zeros(dim)
    U90[0] = -U[1]
    U90[1] = U[0]
    V90[0] = -V[1]
    V90[1] = V[0]

    alpha1 = np.inner(dP, V90) / np.inner(U, V90)
    alpha2 = -np.inner(dP, U90) / np.inner(V, U90)
    P = P1 + alpha1 * U

    return (alpha1, alpha2, P) if with_ratio else P


def normalize(vector):
    tmp = np.sqrt(sum(map(lambda a: a ** 2, vector)))
    if equal(tmp, 0.0):
        return vector

    factor = 1.0 / tmp
    va = np.empty(len(vector), float)
    for i in range(len(vector)):
        va[i] = factor * vector[i]
    return va
