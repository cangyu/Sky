import numpy as np
import math
from scipy.misc import comb
from src.nurbs.basis import *
from src.iges.iges_entity126 import IGES_Entity126


class Curve(object):
    def __init__(self, U, Pw):
        """
        NURBS曲线
        :param U:节点矢量 
        :param Pw: 带权控制点
        """

        self.n = len(Pw) - 1
        self.m = len(U) - 1
        self.p = self.m - self.n - 1

        self.U = np.copy(U)
        self.Pw = np.copy(Pw)

        self.N = Basis(self.U, self.p)

    def w(self, u, d=0):
        ans = 0.0
        for i in range(0, self.n + 1):
            ans += self.N(i, self.p, u, d) * self.Pw[3]

        return ans

    def __call__(self, u, d=0):
        """
        计算曲线上的点
        :param u: 目标参数
        :param d: 求导次数
        :return: 曲线在u处的d阶导矢
        """

        Ad = np.zeros(3)
        nipd = np.zeros(self.n + 1)
        for i in range(0, self.n + 1):
            nipd[i] = self.N(i, p, u, d)
        for i in range(0, self.n + 1):
            for j in range(0, 3):
                Ad[j] += nipd[i] * self.Pw[j]

        ccw = np.zero(3)
        i = 1
        while i <= d:
            ccw += comb(d, i, exact=True) * self.w(u, i) * self.__call__(u, d - i)

        wu = self.w(u)
        if equal(wu, 0.0):
            raise ValueError("Invalid weight settings!")

        return (Ad - ccw) / wu

    def to_iges(self, planar, closed, polynomial, periodic, sp, ep, norm, form=0):

        w = np.zeros(self.n + 1)
        cpt = np.zeros((self.n + 1, 3))

        for i in range(0, self.n + 1):
            w[i] = self.Pw[i][3]
            for j in range(0, 3):
                cpt[i][j] = Pw[i][j] / w[i]

        return IGES_Entity126(self.p, self.n, planar, closed, polynomial, periodic, self.U, w, cpt, sp, ep, norm, form)
