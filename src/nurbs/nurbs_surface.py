import numpy as np
import math
from scipy.misc import comb
from src.nurbs.basis import *
from src.iges.iges_entity128 import IGES_Entity128


class NURBS_Surface(object):
    def __init__(self, U, V, Pw):
        self.U = np.copy(U)
        self.V = np.copy(V)
        self.Pw = np.copy(Pw)

        self.n, self.m, dim = Pw.shape
        self.n = self.n - 1
        self.m = self.m - 1

        self.p = len(U) - 1 - self.n - 1
        self.q = len(V) - 1 - self.m - 1

        self.N1 = Basis(U, self.p)
        self.N2 = Basis(V, self.q)

        self.weight = np.zeros((self.n + 1, self.m + 1))
        self.cpt = np.zeros((self.n + 1, self.m + 1, 3))
        self.isPoly = True

        for i in range(0, self.n + 1):
            for j in range(0, self.m + 1):
                self.weight[i][j] = self.Pw[i][j][3]
                if self.isPoly and (not equal(self.weight[i][j], 1.0)):
                    self.isPoly = False
                for k in range(0, 3):
                    self.cpt[i][j][k] = self.Pw[i][j][k] / self.weight[i][j]

    def w(self, u, v, k=0, l=0):
        ans = 0.0
        nipd = np.zeros(self.n + 1)
        mjqd = np.zeros(self.m + 1)
        for i in range(0, self.n + 1):
            nipd[i] = self.N1(i, p, u, k)
        for j in range(0, self.m + 1):
            mjqd[j] = self.N2(j, q, v, l)

        for i in range(0, self.n + 1):
            for j in range(0, self.m + 1):
                ans += nipd[i] * mjqd[j] * self.Pw[i][j][3]

        return ans

    def __call__(self, u, v, k=0, l=0):
        nipd = np.zeros(self.n + 1)
        mjqd = np.zeros(self.m + 1)
        for i in range(0, self.n + 1):
            nipd[i] = self.N1(i, p, u, k)
        for j in range(0, self.m + 1):
            mjqd[j] = self.N2(j, q, v, l)

        Akl = np.zeros(3)
        for i in range(0, self.n + 1):
            for j in range(0, self.m + 1):
                for kk in range(0, 3):
                    Akl[kk] += nipd[i] * mjqd[j] * self.Pw[i][j][kk]

        ccw1 = np.zeros(3)
        i = 1
        while i <= k:
            ccw1 += comb(k, i, exact=True) * self.w(u, v, i, 0) * self.__call__(u, v, k - i, l)
            i += 1

        ccw2 = np.zeros(3)
        j = 1
        while j <= l:
            ccw2 += comb(l, j, exact=True) * self.w(u, v, 0, j) * self.__call__(u, v, k, l - j)
            j += 1

        ccw3 = np.zeros(3)
        i = 1
        while i <= k:
            ci = comb(k, i, exact=True)
            j = 1
            while j <= l:
                ccw3 += ci * comb(l, j, exact=True) * self.w(u, v, i, j) * self.__call__(u, v, k - i, l - j)
                j += 1
            i += 1

        wuv = self.w(u, v)
        if equal(wuv, 0.0):
            raise ValueError("Invalid weight settings!")

        return (Akl - ccw1 - ccw2 - ccw3) / wuv

    def to_iges(self, closed_u, closed_v, periodic_u, periodic_v, form=0):
        return IGES_Entity128(self.U, self.V, self.p, self.q, self.n, self.m, self.weight, self.cpt,
                              closed_u, closed_v, (1 if self.isPoly else 0), periodic_u, periodic_v,
                              self.U[0], self.U[-1], self.V[0], self.V[-1], form)
