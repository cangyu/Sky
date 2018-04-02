#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from abc import abstractmethod
import math
import numpy as np
from scipy.linalg import norm
from scipy import sparse
from scipy.sparse.linalg import dsolve
from misc import vector_square

"""
Implementation of the grid smoothing tools using elliptic PDE.

Note:
1. (i,j,k) is corresponding to (x,y,z),(u,v.w),(xi, eta, zeta),(x1,x2,x3),(I,J,K)...
2. (N,M) notation is not suggested to avoid confuse on column-and-row.
3. By default, delta_xi = delta_eta = delta_zeta = 1, and is neglected in code.
4. All the partial derivatives in elliptic PDE are calculated with the central scheme.
"""


class EllipticGrid2D(object):
    di = [0, 1, 0, -1, 0, 1, -1, -1, 1]
    dj = [0, 0, 1, 0, -1, 1, 1, -1, -1]

    def __init__(self, grid):
        """
        2D curvilinear grid based on the elliptic PDE.
        :param grid: Initial grid. The subscript iterate through (Dim1, Dim2, Dim3), each element contains (X, Y, Z).
        """

        '''Pre-check'''
        assert len(grid.shape) == 3
        assert grid.shape[-1] in (2, 3)

        '''Shape constants'''
        ii, jj, dim = grid.shape
        assert (ii - 2) * (jj - 2) != 0

        '''Grid'''
        self.r = np.copy(grid[:, :, :2])

        '''Partial derivatives'''
        self.r1 = np.zeros((ii, jj, 2))
        self.r2 = np.zeros((ii, jj, 2))
        self.r11 = np.zeros((ii, jj, 2))
        self.r22 = np.zeros((ii, jj, 2))
        self.r12 = np.zeros((ii, jj, 2))

        '''Source term'''
        self.pq = np.zeros((ii, jj, 2))

        '''Coefficients'''
        self.a = np.zeros((ii, jj))  # alpha
        self.b = np.zeros((ii, jj))  # beta
        self.g = np.zeros((ii, jj))  # gamma
        self.j2 = np.zeros((ii, jj))  # square of det(jacobi)

    @property
    def i_num(self):
        return self.r.shape[0]

    @property
    def j_num(self):
        return self.r.shape[1]

    def is_special(self, i, j):
        return i == 0 or i == self.i_num - 1 or j == 0 or j == self.j_num - 1

    def internal_pnt_idx(self, i, j):
        return (i - 1) + (j - 1) * (self.i_num - 2)

    @property
    def grid(self):
        return self.r

    def r_xi(self, i, j):
        return 0.5 * (self.r[i + 1][j] - self.r[i - 1][j])

    def r_eta(self, i, j):
        return 0.5 * (self.r[i][j + 1] - self.r[i][j - 1])

    def r_xi2(self, i, j):
        return self.r[i + 1][j] - 2 * self.r[i][j] + self.r[i - 1][j]

    def r_eta2(self, i, j):
        return self.r[i][j + 1] - 2 * self.r[i][j] + self.r[i][j - 1]

    def r_xi_eta(self, i, j):
        return 0.25 * (self.r[i + 1][j + 1] + self.r[i - 1][j - 1] - self.r[i + 1][j - 1] - self.r[i - 1][j + 1])

    def alpha(self, i, j):
        return vector_square(self.r_eta(i, j))

    def beta(self, i, j):
        return np.dot(self.r_eta(i, j), self.r_xi(i, j))

    def gamma(self, i, j):
        return vector_square(self.r_xi(i, j))

    def jacobi(self, i, j):
        return np.linalg.det(np.matrix([self.r_xi(i, j), self.r_eta(i, j)])) ** 2

    @abstractmethod
    def calc_all_param(self):
        pass

    @abstractmethod
    def calc_eqn_param(self, i, j):
        pass

    def smooth(self):
        """
        Smooth the grid with Picard iteration.
        :return: None.
        """

        var_num = (self.i_num - 2) * (self.j_num - 2)
        rhs = np.empty((var_num, 2))

        '''Solve the grid iteratively'''
        iteration_cnt = 0
        residual = sys.float_info.max
        while not math.isclose(residual, 0, abs_tol=1e-5):
            '''Calculate all coefficients and derivatives'''
            self.calc_all_param()

            '''Build Ax = b'''
            rhs.fill(0.0)
            eqn_idx = 0
            row = []
            col = []
            val = []
            for j in range(1, self.j_num - 1):
                for i in range(1, self.i_num - 1):
                    ca = self.calc_eqn_param(i, j)  # surrounding coefficients
                    for t in range(9):
                        ii = i + EllipticGrid2D.di[t]
                        jj = j + EllipticGrid2D.dj[t]
                        if self.is_special(ii, jj):
                            rhs[eqn_idx] -= ca[t] * self.r[ii][jj]
                        else:
                            row.append(eqn_idx)
                            col.append(self.internal_pnt_idx(ii, jj))
                            val.append(ca[t])
                    eqn_idx += 1

            '''Construct the sparse coefficient matrix'''
            scm = sparse.coo_matrix((val, (row, col)), shape=(var_num, var_num), dtype=float).tocsr()

            '''Solve the grid'''
            u = np.copy([dsolve.spsolve(scm, b) for b in rhs.transpose()]).transpose()

            '''Update and calculate residual'''
            eqn_idx = 0
            residual = sys.float_info.min
            for j in range(1, self.j_num - 1):
                for i in range(1, self.i_num - 1):
                    cur_residual = norm(u[eqn_idx] - self.r[i][j], np.inf)
                    if cur_residual > residual:
                        residual = cur_residual
                    self.r[i][j] = u[eqn_idx]
                    eqn_idx += 1

            iteration_cnt += 1


class Laplace2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        Smooth the grid with Laplace PDE.
        :param grid: The initial grid.
        """

        super(Laplace2D, self).__init__(grid)

    def calc_all_param(self):
        for i in range(1, self.i_num - 1):
            for j in range(1, self.j_num - 1):
                self.r1[i][j] = self.r_xi(i, j)
                self.r2[i][j] = self.r_eta(i, j)
                self.a[i][j] = vector_square(self.r2[i][j])
                self.b[i][j] = np.dot(self.r1[i][j], self.r2[i][j])
                self.g[i][j] = vector_square(self.r1[i][j])

    def calc_eqn_param(self, i, j):
        ans = np.empty(9, float)
        ans[0] = -2 * (self.a[i][j] + self.g[i][j])
        ans[1] = ans[3] = self.a[i][j]
        ans[2] = ans[4] = self.g[i][j]
        ans[6] = ans[8] = self.b[i][j] / 2
        ans[5] = ans[7] = -ans[6]
        return ans


class ThomasMiddlecoff2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        Smooth the grid with the Thomas-Middlecoff method.
        :param grid: Initial grid.
        """

        super(ThomasMiddlecoff2D, self).__init__(grid)

        self.phi = np.zeros((self.i_num, self.j_num))
        self.psi = np.zeros((self.i_num, self.j_num))

    def calc_all_param(self):
        for i in range(1, self.i_num - 1):
            for j in range(1, self.j_num - 1):
                self.r1[i][j] = self.r_xi(i, j)
                self.r2[i][j] = self.r_eta(i, j)
                self.r11[i][j] = self.r_xi2(i, j)
                self.r22[i][j] = self.r_eta2(i, j)
                self.r12[i][j] = self.r_xi_eta(i, j)
                self.a[i][j] = vector_square(self.r2[i][j])
                self.b[i][j] = np.dot(self.r1[i][j], self.r2[i][j])
                self.g[i][j] = vector_square(self.r1[i][j])
        for j in (0, self.j_num - 1):
            for i in range(1, self.i_num - 1):
                self.r1[i][j] = self.r_xi(i, j)
                self.r11[i][j] = self.r_xi2(i, j)
                self.phi[i][j] = - np.dot(self.r1[i][j], self.r11[i][j]) / vector_square(self.r1[i][j])
        for i in (0, self.i_num - 1):
            for j in range(1, self.j_num - 1):
                self.r2[i][j] = self.r_eta(i, j)
                self.r22[i][j] = self.r_eta2(i, j)
                self.psi[i][j] = -np.dot(self.r2[i][j], self.r22[i][j]) / vector_square(self.r2[i][j])
        for i in range(1, self.i_num - 1):
            dist = np.linspace(self.phi[i][0], self.phi[i][self.j_num - 1], self.j_num)
            for j in range(1, self.j_num - 1):
                self.phi[i][j] = dist[j]
        for j in range(1, self.j_num - 1):
            dist = np.linspace(self.psi[0][j], self.psi[self.i_num - 1][j], self.i_num)
            for i in range(1, self.i_num - 1):
                self.psi[i][j] = dist[i]

    def calc_eqn_param(self, i, j):
        ans = np.empty(9, float)
        ans[0] = -2.0 * (self.a[i][j] + self.g[i][j])
        ans[1] = self.a[i][j] * (1 + self.phi[i][j] / 2)
        ans[2] = self.g[i][j] * (1 + self.psi[i][j] / 2)
        ans[3] = self.a[i][j] * (1 - self.phi[i][j] / 2)
        ans[4] = self.g[i][j] * (1 - self.psi[i][j] / 2)
        ans[6] = ans[8] = self.b[i][j] / 2
        ans[5] = ans[7] = -ans[6]
        return ans


class EllipticGrid3D(object):
    di = [0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1]
    dj = [0, 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1, 0, 0, 0, 0]
    dk = [0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1]

    def __init__(self, grid):
        """
        Base class for 3D elliptic-PDE based grid smoother.
        :param grid: Initial grid.
        """

        '''Pre-check'''
        assert len(grid.shape) == 4
        assert grid.shape[-1] == 3

        '''Shape constants'''
        ii, jj, kk, dim = grid.shape
        assert (ii - 2) * (jj - 2) * (kk - 2) != 0

        '''Grid'''
        self.r = np.copy(grid)

        '''Derivatives'''
        self.r1 = np.zeros_like(self.r)
        self.r2 = np.zeros_like(self.r)
        self.r3 = np.zeros_like(self.r)
        self.r11 = np.zeros_like(self.r)
        self.r22 = np.zeros_like(self.r)
        self.r33 = np.zeros_like(self.r)
        self.r12 = np.zeros_like(self.r)
        self.r23 = np.zeros_like(self.r)
        self.r31 = np.zeros_like(self.r)

        '''Source term'''
        self.pqr = np.zeros((ii, jj, kk, 3))

        '''Equation coefficients'''
        self.a1 = np.zeros((ii, jj, kk))  # alpha1
        self.a2 = np.zeros((ii, jj, kk))  # alpha2
        self.a3 = np.zeros((ii, jj, kk))  # alpha3
        self.b12 = np.zeros((ii, jj, kk))  # beta12
        self.b23 = np.zeros((ii, jj, kk))  # beta23
        self.b31 = np.zeros((ii, jj, kk))  # beta23
        self.j2 = np.zeros((ii, jj, kk))  # Jacobi

    @property
    def i_num(self):
        return self.r.shape[0]

    @property
    def j_num(self):
        return self.r.shape[1]

    @property
    def k_num(self):
        return self.r.shape[2]

    def is_special(self, i, j, k):
        return i == 0 or i == self.i_num - 1 or j == 0 or j == self.j_num - 1 or k == 0 or k == self.k_num - 1

    def internal_pnt_idx(self, i, j, k):
        return (k - 1) * (self.j_num - 2) * (self.i_num - 2) + (j - 1) * (self.i_num - 2) + (i - 1)

    @property
    def grid(self):
        return self.r

    def r_xi(self, i, j, k):
        return 0.5 * (self.r[i + 1][j][k] - self.r[i - 1][j][k])

    def r_eta(self, i, j, k):
        return 0.5 * (self.r[i][j + 1][k] - self.r[i][j - 1][k])

    def r_zeta(self, i, j, k):
        return 0.5 * (self.r[i][j][k + 1] - self.r[i][j][k - 1])

    def r_xi2(self, i, j, k):
        return self.r[i - 1][j][k] - 2 * self.r[i][j][k] + self.r[i + 1][j][k]

    def r_eta2(self, i, j, k):
        return self.r[i][j - 1][k] - 2 * self.r[i][j][k] + self.r[i][j + 1][k]

    def r_zeta2(self, i, j, k):
        return self.r[i][j][k - 1] - 2 * self.r[i][j][k] + self.r[i][j][k + 1]

    def r_xi_eta(self, i, j, k):
        return 0.25 * (self.r[i + 1][j + 1][k] - self.r[i + 1][j - 1][k] - self.r[i - 1][j + 1][k] + self.r[i - 1][j - 1][k])

    def r_eta_zeta(self, i, j, k):
        return 0.25 * (self.r[i][j + 1][k + 1] - self.r[i][j - 1][k + 1] - self.r[i][j + 1][k - 1] + self.r[i][j - 1][k - 1])

    def r_zeta_xi(self, i, j, k):
        return 0.25 * (self.r[i + 1][j][k + 1] - self.r[i - 1][j][k + 1] - self.r[i + 1][j][k - 1] + self.r[i - 1][j][k - 1])

    def alpha1(self, i, j, k):
        r2 = self.r_eta(i, j, k)
        r3 = self.r_zeta(i, j, k)
        return vector_square(r2) * vector_square(r3) - np.dot(r2, r3) ** 2

    def alpha2(self, i, j, k):
        r3 = self.r_zeta(i, j, k)
        r1 = self.r_xi(i, j, k)
        return vector_square(r3) * vector_square(r1) - np.dot(r3, r1) ** 2

    def alpha3(self, i, j, k):
        r1 = self.r_xi(i, j, k)
        r2 = self.r_eta(i, j, k)
        return vector_square(r1) * vector_square(r2) - np.dot(r1, r2) ** 2

    def beta12(self, i, j, k):
        r1 = self.r_xi(i, j, k)
        r2 = self.r_eta(i, j, k)
        r3 = self.r_zeta(i, j, k)
        return np.dot(r1, r3) * np.dot(r2, r3) - np.dot(r1, r2) * vector_square(r3)

    def beta23(self, i, j, k):
        r1 = self.r_xi(i, j, k)
        r2 = self.r_eta(i, j, k)
        r3 = self.r_zeta(i, j, k)
        return np.dot(r2, r1) * np.dot(r3, r1) - np.dot(r2, r3) * vector_square(r1)

    def beta31(self, i, j, k):
        r1 = self.r_xi(i, j, k)
        r2 = self.r_eta(i, j, k)
        r3 = self.r_zeta(i, j, k)
        return np.dot(r3, r2) * np.dot(r1, r2) - np.dot(r3, r1) * vector_square(r2)

    def jacobi(self, i, j, k):
        return np.linalg.det(np.matrix([self.r_xi(i, j, k), self.r_eta(i, j, k), self.r_zeta(i, j, k)])) ** 2

    @abstractmethod
    def calc_all_param(self):
        pass

    @abstractmethod
    def calc_eqn_param(self, i, j, k):
        pass

    def smooth(self):
        """
        Smooth the grid using Picard Iteration.
        :return: None.
        """

        var_num = (self.i_num - 2) * (self.j_num - 2) * (self.k_num - 2)
        rhs = np.zeros((var_num, 3))

        '''Solve the grid iteratively'''
        iteration_cnt = 0
        residual = sys.float_info.max
        while not math.isclose(residual, 0, abs_tol=1e-5):
            '''Calculate all coefficients'''
            self.calc_all_param()

            '''Build Ax=b'''
            eqn_idx = 0
            rhs.fill(0.0)
            row = []
            col = []
            val = []
            for k in range(1, self.k_num - 1):
                for j in range(1, self.j_num - 1):
                    for i in range(1, self.i_num - 1):
                        '''Calculate stencil coefficients'''
                        ca = self.calc_eqn_param(i, j, k)

                        '''Construct the equation'''
                        for t in range(19):
                            ii = i + EllipticGrid3D.di[t]
                            jj = j + EllipticGrid3D.dj[t]
                            kk = k + EllipticGrid3D.dk[t]
                            if self.is_special(ii, jj, kk):
                                rhs[eqn_idx] -= ca[t] * self.r[ii][jj][kk]
                            else:
                                row.append(eqn_idx)
                                col.append(self.internal_pnt_idx(ii, jj, kk))
                                val.append(ca[t])

                        eqn_idx += 1

            '''Construct the sparse coefficient matrix'''
            scm = sparse.coo_matrix((val, (row, col)), shape=(var_num, var_num), dtype=float).tocsr()

            '''Solve the grid'''
            u = np.copy([dsolve.spsolve(scm, b) for b in rhs.transpose()]).transpose()

            '''Update and calculate residual'''
            eqn_idx = 0
            residual = sys.float_info.min
            for k in range(1, self.k_num - 1):
                for j in range(1, self.j_num - 1):
                    for i in range(1, self.i_num - 1):
                        cur_residual = norm(u[eqn_idx] - self.r[i][j][k], np.inf)
                        if cur_residual > residual:
                            residual = cur_residual
                        self.r[i][j][k] = u[eqn_idx]
                        eqn_idx += 1

            iteration_cnt += 1
            print(residual)


class Laplace3D(EllipticGrid3D):
    def __init__(self, grid):
        super(Laplace3D, self).__init__(grid)

    def calc_all_param(self):
        for i in range(1, self.i_num - 1):
            for j in range(1, self.j_num - 1):
                for k in range(1, self.k_num - 1):
                    self.r1[i][j][k] = self.r_xi(i, j, k)
                    self.r2[i][j][k] = self.r_eta(i, j, k)
                    self.r3[i][j][k] = self.r_zeta(i, j, k)
                    self.r11[i][j][k] = self.r_xi2(i, j, k)
                    self.r22[i][j][k] = self.r_eta2(i, j, k)
                    self.r33[i][j][k] = self.r_zeta2(i, j, k)
                    self.r12[i][j][k] = self.r_xi_eta(i, j, k)
                    self.r23[i][j][k] = self.r_eta_zeta(i, j, k)
                    self.r31[i][j][k] = self.r_zeta_xi(i, j, k)
                    self.a1[i][j][k] = vector_square(self.r2[i][j][k]) * vector_square(self.r3[i][j][k]) - np.dot(self.r2[i][j][k], self.r3[i][j][k]) ** 2
                    self.a2[i][j][k] = vector_square(self.r3[i][j][k]) * vector_square(self.r1[i][j][k]) - np.dot(self.r3[i][j][k], self.r1[i][j][k]) ** 2
                    self.a3[i][j][k] = vector_square(self.r1[i][j][k]) * vector_square(self.r2[i][j][k]) - np.dot(self.r1[i][j][k], self.r2[i][j][k]) ** 2
                    self.b12[i][j][k] = np.dot(self.r1[i][j][k], self.r3[i][j][k]) * np.dot(self.r2[i][j][k], self.r3[i][j][k]) - np.dot(self.r1[i][j][k], self.r2[i][j][k]) * vector_square(self.r3[i][j][k])
                    self.b23[i][j][k] = np.dot(self.r2[i][j][k], self.r1[i][j][k]) * np.dot(self.r3[i][j][k], self.r1[i][j][k]) - np.dot(self.r2[i][j][k], self.r3[i][j][k]) * vector_square(self.r1[i][j][k])
                    self.b31[i][j][k] = np.dot(self.r3[i][j][k], self.r2[i][j][k]) * np.dot(self.r1[i][j][k], self.r2[i][j][k]) - np.dot(self.r3[i][j][k], self.r1[i][j][k]) * vector_square(self.r2[i][j][k])

    def calc_eqn_param(self, i, j, k):
        ans = np.empty(19, float)
        ans[0] = -2 * (self.a1[i][j][k] + self.a2[i][j][k] + self.a3[i][j][k])
        ans[1] = ans[2] = self.a1[i][j][k]
        ans[3] = ans[4] = self.a2[i][j][k]
        ans[5] = ans[6] = self.a3[i][j][k]
        ans[7] = ans[8] = 0.5 * self.b12[i][j][k]
        ans[9] = ans[10] = -ans[8]
        ans[11] = ans[12] = 0.5 * self.b23[i][j][k]
        ans[13] = ans[14] = -ans[12]
        ans[15] = ans[16] = 0.5 * self.b31[i][j][k]
        ans[17] = ans[18] = -ans[16]
        return ans


class ThomasMiddlecoff3D(EllipticGrid3D):
    def __init__(self, grid):
        super(ThomasMiddlecoff3D, self).__init__(grid)

        self.phi = np.zeros((self.i_num, self.j_num, self.k_num))
        self.psi = np.zeros((self.i_num, self.j_num, self.k_num))
        self.omega = np.zeros((self.i_num, self.j_num, self.k_num))

    def calc_all_param(self):
        pass

    def calc_eqn_param(self, i, j, k):
        pass
