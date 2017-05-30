import numpy as np
import math


def equal(a: float, b: float):
    """
    判断两个浮点数是否相等
    """
    return (True if math.fabs(a - b) < 1e-8 else False)


def PntDist(lhs, rhs):
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

    pw = np.ones(len(pnt) + 1) * w
    for i in range(0, len(pnt)):
        pw[i] = pnt[i] if equal(w, 0.0) else pnt[i] * w

    return pw


def to_cartesian(pnt):
    """
    转换为笛卡尔坐标
    :param pnt: 原始齐次坐标
    :return: 投影后得到的笛卡尔坐标
    """
    p = np.zeros(len(pnt) - 1)
    for i in range(0, len(p)):
        p[i] = pnt[i] if equal(pnt[-1], 0.0) else pnt[i] / pnt[-1]

    return p


def find_span(U, u):
    mi = U[0]
    ma = U[-1]
    left = 0
    right = len(U) - 1
    while U[left + 1] == mi:
        left += 1
    while U[right - 1] == ma:
        right -= 1

    for i in range(left, right):
        if U[i] <= u < U[i + 1]:
            return i

    return right - 1  # Corner case when u=U[-1]


def all_basis_funs(i, u, p, U):
    N = np.zeros(p + 1)
    left = np.zeros(p + 1)
    right = np.zeros(p + 1)
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
    ans = np.zeros(n + 1)
    for k in range(i - p, i + 1):
        ans[k] = N[k - (i - p)]

    return ans

class Basis(object):
    def __init__(self, U, p):
        """
        定义在非周期节点矢量上的B样条基函数
        :param U: 节点矢量
        :param p: 基函数次数
        """

        self.U = np.copy(U)
        self.p = p
        self.m = len(self.U) - 1
        self.n = self.m - self.p - 1

    def __call__(self, i, p, u, d=0):
        """
        计算基函数导数
        :param i: 第i个B样条基函数
        :param p: B样条基函数次数
        :param u: 目标参数
        :param d: 求导次数
        :return: N(i, p, u)的d阶导数
        """

        if d < 0:
            raise ValueError("Order of Derivative should be non-negative!")

        if d == 0:
            if p == 0:
                if equal(u, self.U[-1]):
                    return 1.0 if i == self.n else 0.0
                else:
                    return 1.0 if self.U[i] <= u < self.U[i + 1] else 0.0
            else:
                ans = 0.0
                tmp = self.U[i + p] - self.U[i]
                if not equal(tmp, 0.0):
                    ans += (u - self.U[i]) / tmp * self.__call__(i, p - 1, u)

                tmp = self.U[i + p + 1] - self.U[i + 1]
                if not equal(tmp, 0.0):
                    ans += (self.U[i + p + 1] - u) / tmp * self.__call__(i + 1, p - 1, u)

                return ans
        else:
            ans = 0.0
            tmp = self.U[i + p] - self.U[i]
            if not equal(tmp, 0.0):
                ans += self.__call__(i, p - 1, u, d - 1) / tmp

            tmp = self.U[i + p + 1] - self.U[i + 1]
            if not equal(tmp, 0.0):
                ans += self.__call__(i + 1, p - 1, u, d - 1) / tmp

            ans *= p

            return ans
