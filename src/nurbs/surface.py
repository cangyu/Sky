import numpy as np
import math
from src.nurbs.basis import *


class Surface(object):
    def __init__(self, U, V, Pw):
        self.U = np.copy(U)
        self.V = np.copy(V)
        self.Pw = np.copy(Pw)

        self.n, self.m, dim = Pw.shape
        self.n = self.n - 1
        self.m = self.m - 1

        self.p = self.len(U) - 1 - self.n - 1
        self.q = self.len(V) - 1 - self.m - 1

        self.N1 = Basis(U, self.p)
        self.N2 = Basis(V, self.q)

    def __call__(self, u, v, k, l):
        pass
