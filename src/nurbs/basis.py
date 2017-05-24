import numpy as np
import math


def equal(a: float, b: float):
    return True if math.fabs(a - b) < 1e-8 else False


class Basis(object):
    def __init__(self, U, p):
        """
        定义在非周期节点矢量上的B样条基函数
        :param U: 节点矢量
        :param p: 基函数次数
        """

        self.U = np.copy(U)

        self.p = p
        self.m = len(self.U)
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
            raise ValueError("Order of Derivative shold be non-negative!")

        if d == 0:
            if p == 0:
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
