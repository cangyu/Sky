import unittest
import sys
import math
import numpy as np
from abc import abstractmethod
from scipy.linalg import norm
from scipy import sparse
from scipy.sparse.linalg import dsolve
from misc import vector_square

"""
Implementation of the grid smoothing tools using elliptic PDE operator.

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
        :param grid: Initial grid. The subscript iterate through (Dim1, Dim2), each element contains (X, Y).
        """

        '''Pre-check'''
        if len(grid.shape) != 3:
            raise AssertionError("Invalid input grid.")
        if grid.shape[-1] != 2:
            raise AssertionError('Invalid dimension.')

        '''Shape constants'''
        ii, jj, dim = grid.shape

        '''Grid'''
        self.r = np.copy(grid)

        '''Partial derivatives'''
        self.r1 = np.zeros_like(grid)
        self.r2 = np.zeros_like(grid)
        self.r11 = np.zeros_like(grid)
        self.r22 = np.zeros_like(grid)
        self.r12 = np.zeros_like(grid)

        '''Source term'''
        self.pq = np.zeros_like(grid)

        '''Coefficients'''
        self.a = np.zeros((ii, jj))  # alpha
        self.b = np.zeros((ii, jj))  # beta
        self.g = np.zeros((ii, jj))  # gamma
        self.j = np.zeros((ii, jj))  # jacobi

    @property
    def i_dim(self):
        return self.r.shape[0]

    @property
    def j_dim(self):
        return self.r.shape[1]

    def is_special(self, i, j):
        return i == 0 or i == self.i_dim - 1 or j == 0 or j == self.j_dim - 1

    def internal_pnt_idx(self, i, j):
        return (i - 1) + (j - 1) * (self.i_dim - 2)

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
        return vector_square(self.r2[i][j])

    def beta(self, i, j):
        return np.dot(self.r2[i][j], self.r1[i][j])

    def gamma(self, i, j):
        return vector_square(self.r1[i][j])

    def jacobi(self, i, j):
        return np.linalg.det(np.matrix([self.r1[i][j], self.r2[i][j]]))

    @abstractmethod
    def smooth(self):
        pass


class Laplace2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        Smooth the grid with Laplace operator.
        :param grid: The initial grid.
        """

        super(Laplace2D, self).__init__(grid)

    def calc_all_param(self):
        for i in range(1, self.i_dim - 1):
            for j in range(1, self.j_dim - 1):
                self.r1[i][j] = self.r_xi(i, j)
                self.r2[i][j] = self.r_eta(i, j)
                self.r11[i][j] = self.r_xi2(i, j)
                self.r22[i][j] = self.r_eta2(i, j)
                self.r12[i][j] = self.r_xi_eta(i, j)
                self.a[i][j] = self.alpha(i, j)
                self.b[i][j] = self.beta(i, j)
                self.g[i][j] = self.gamma(i, j)

    def calc_eqn_param(self, i, j):
        ans = np.empty(9)
        ans[0] = -2 * (self.a[i][j] + self.g[i][j])
        ans[1] = ans[3] = self.a[i][j]
        ans[2] = ans[4] = self.g[i][j]
        ans[5] = ans[7] = -self.b[i][j] / 2
        ans[6] = ans[8] = -ans[5]
        return ans

    def smooth(self):
        """
        Smooth the grid with Picard iteration.
        :return: None.
        """

        '''Temp vars used to construct the sparse coefficient matrix'''
        var_num = (self.i_dim - 1) * (self.j_dim - 1)
        rhs = np.empty((var_num, 2))

        '''Solve the grid iteratively'''
        iteration_cnt = 0
        residual = sys.float_info.max
        while not math.isclose(residual, 0, abs_tol=1e-5):
            '''Calculate all coefficients and derivatives'''
            self.calc_all_param()

            '''Build Ar = b'''
            rhs.fill(0.0)
            eqn_idx = 0
            row = []
            col = []
            val = []
            for i in range(1, self.i_dim - 1):
                for j in range(1, self.j_dim - 1):
                    ca = self.calc_eqn_param(i, j)  # surrounding coefficients

                    '''Construct the equation'''
                    for t in range(9):
                        ii = i + EllipticGrid2D.di[t]
                        jj = j + EllipticGrid2D.dj[t]
                        if self.is_special(i, j):
                            rhs[eqn_idx] -= ca[t] * self.r[ii][jj]
                        else:
                            row.append(eqn_idx)
                            col.append(self.internal_pnt_idx(ii, jj))
                            val.append(ca[t])

                    eqn_idx += 1

            '''Construct the sparse coefficient matrix'''
            scm = sparse.coo_matrix((val, (row, col)), shape=(var_num, var_num), dtype=float).tocsr()

            '''Solve the grid'''
            u = dsolve.spsolve(scm, rhs)
            eqn_idx = 0
            for i in range(1, self.i_dim - 1):
                for j in range(1, self.j_dim - 1):
                    residual = max(residual, norm(u[eqn_idx] - self.r[i][j], np.inf))
                    eqn_idx += 1

            iteration_cnt += 1


class ThomasMiddlecoff2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        Smooth the grid with the Thomas-Middlecoff method.
        :param grid: Initial grid.
        """

        super(ThomasMiddlecoff2D, self).__init__(grid)

        self.phi = np.zeros((self.i_dim, self.j_dim))
        self.psi = np.zeros((self.i_dim, self.j_dim))

    def calc_all_param(self):
        for i in range(1, self.i_dim - 1):
            for j in range(1, self.j_dim - 1):
                self.r1[i][j] = self.r_xi(i, j)
                self.r2[i][j] = self.r_eta(i, j)
                self.r11[i][j] = self.r_xi2(i, j)
                self.r22[i][j] = self.r_eta2(i, j)
                self.r12[i][j] = self.r_xi_eta(i, j)
                self.a[i][j] = self.alpha(i, j)
                self.b[i][j] = self.beta(i, j)
                self.g[i][j] = self.gamma(i, j)

        '''phi, psi: Boundary first, then interpolate internal space linearly'''
        for j in (0, self.j_dim - 1):
            for i in range(1, self.i_dim - 1):
                self.phi[i][j] = - np.dot(self.r1[i][j], self.r11[i][j]) / vector_square(self.r1[i][j])
        for i in (0, self.i_dim - 1):
            for j in range(1, self.j_dim - 1):
                self.psi[i][j] = -np.dot(self.r2[i][j], self.r22[i][j]) / vector_square(self.r2[i][j])
        for i in range(1, self.i_dim - 1):
            dist = np.linspace(self.phi[i][0], self.phi[i][self.j_dim - 1], self.j_dim)
            for j in range(1, self.j_dim - 1):
                self.phi[i][j] = dist[j]
        for j in range(1, self.j_dim - 1):
            dist = np.linspace(self.psi[0][j], self.psi[self.i_dim - 1][j], self.i_dim)
            for i in range(1, self.i_dim - 1):
                self.psi[i][j] = dist[i]

    def calc_eqn_param(self, i, j):
        ans = np.empty(9, float)
        ans[0] = -2.0 * (self.a[i][j] + self.g[i][j])
        ans[1] = self.a[i][j] * (1 + self.phi[i][j] / 2)
        ans[2] = self.g[i][j] * (1 + self.psi[i][j] / 2)
        ans[3] = self.a[i][j] * (1 - self.phi[i][j] / 2)
        ans[4] = self.g[i][j] * (1 - self.psi[i][j] / 2)
        ans[5] = ans[7] = -self.b[i][j] / 2
        ans[6] = ans[8] = -ans[5]
        return ans

    def smooth(self):
        var_num = (self.i_dim - 1) * (self.j_dim - 1)
        rhs = np.empty((var_num, 2))
        iteration_cnt = 0
        residual = sys.float_info.max
        while not math.isclose(residual, 0, abs_tol=1e-5):
            '''Calculate all coefficients and derivatives'''
            self.calc_all_param()

            '''Build Ar = b'''
            eqn_idx = 0
            row = []
            col = []
            val = []
            rhs.fill(0.0)
            for i in range(1, self.i_dim - 1):
                for j in range(1, self.j_dim - 1):
                    ca = self.calc_eqn_param(i, j)
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
            u = dsolve.spsolve(scm, rhs)
            eqn_idx = 0
            for i in range(1, self.i_dim - 1):
                for j in range(1, self.j_dim - 1):
                    residual = max(residual, norm(u[eqn_idx] - self.r[i][j], np.inf))
                    eqn_idx += 1

            iteration_cnt += 1


class EllipticGrid3D(object):
    di = [0, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1]
    dj = [0, 0, 0, -1, 1, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 0, 0, 0, 0]
    dk = [0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, -1, 1]

    def __init__(self, grid):
        """
        Base class for 3D elliptic-PDE based grid smoother.
        :param grid: Initial grid.
        """

        '''Pre-check'''
        if len(grid.shape) != 4 or grid.shape[-1] > 3:
            raise AssertionError("Invalid input grid!")

        '''Shape constants'''
        ii, jj, kk, dim = grid.shape

        '''Grid'''
        self.r = np.copy(grid)

        '''Derivatives'''
        self.r1 = np.zeros((ii, jj, kk, 3))
        self.r2 = np.zeros((ii, jj, kk, 3))
        self.r3 = np.zeros((ii, jj, kk, 3))
        self.r11 = np.zeros((ii, jj, kk, 3))
        self.r22 = np.zeros((ii, jj, kk, 3))
        self.r33 = np.zeros((ii, jj, kk, 3))
        self.r12 = np.zeros((ii, jj, kk, 3))
        self.r23 = np.zeros((ii, jj, kk, 3))
        self.r31 = np.zeros((ii, jj, kk, 3))

        '''Source term'''
        self.pqr = np.zeros((ii, jj, kk, 3))

        '''Equation coefficients'''
        self.a1 = np.zeros((ii, jj, kk))  # alpha1
        self.a2 = np.zeros((ii, jj, kk))  # alpha2
        self.a3 = np.zeros((ii, jj, kk))  # alpha3
        self.b12 = np.zeros((ii, jj, kk))  # beta12
        self.b23 = np.zeros((ii, jj, kk))  # beta23
        self.b31 = np.zeros((ii, jj, kk))  # beta23
        self.j = np.zeros((ii, jj, kk))  # Jacobi

    @property
    def i_dim(self):
        return self.r.shape[0]

    @property
    def j_dim(self):
        return self.r.shape[1]

    @property
    def k_dim(self):
        return self.r.shape[2]

    def is_special(self, i, j, k):
        return i == 0 or i == self.i_dim - 1 or j == 0 or j == self.j_dim - 1 or k == 0 or k == self.k_dim - 1

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
        return vector_square(self.r2[i][j][k]) * vector_square(self.r3[i][j][k]) - np.dot(self.r2[i][j][k], self.r3[i][j][k]) ** 2

    def alpha2(self, i, j, k):
        return vector_square(self.r3[i][j][k]) * vector_square(self.r1[i][j][k]) - np.dot(self.r3[i][j][k], self.r1[i][j][k]) ** 2

    def alpha3(self, i, j, k):
        return vector_square(self.r1[i][j][k]) * vector_square(self.r2[i][j][k]) - np.dot(self.r1[i][j][k], self.r2[i][j][k]) ** 2

    def beta12(self, i, j, k):
        return np.dot(self.r1[i][j][k], self.r3[i][j][k]) * np.dot(self.r2[i][j][k], self.r3[i][j][k]) - np.dot(self.r1[i][j][k], self.r2[i][j][k]) * norm(self.r3[i][j][k])

    def beta23(self, i, j, k):
        return np.dot(self.r2[i][j][k], self.r1[i][j][k]) * np.dot(self.r3[i][j][k], self.r1[i][j][k]) - np.dot(self.r2[i][j][k], self.r3[i][j][k]) * norm(self.r1[i][j][k])

    def beta31(self, i, j, k):
        return np.dot(self.r3[i][j][k], self.r2[i][j][k]) * np.dot(self.r1[i][j][k], self.r2[i][j][k]) - np.dot(self.r3[i][j][k], self.r1[i][j][k]) * norm(self.r2[i][j][k])

    def jacobi(self, i, j, k):
        return np.linalg.det(np.matrix([self.r1[i][j][k], self.r2[i][j][k], self.r3[i][j][k]]))

    @abstractmethod
    def smooth(self):
        pass


class Laplace3D(EllipticGrid3D):
    def __init__(self, grid):
        super(Laplace3D, self).__init__(grid)

    def calc_all_param(self):
        for i in range(1, self.i_dim - 1):
            for j in range(1, self.j_dim - 1):
                for k in range(1, self.k_dim - 1):
                    self.r1[i][j][k] = self.r_xi(i, j, k)
                    self.r2[i][j][k] = self.r_eta(i, j, k)
                    self.r3[i][j][k] = self.r_zeta(i, j, k)
                    self.r11[i][j][k] = self.r_xi2(i, j, k)
                    self.r22[i][j][k] = self.r_eta2(i, j, k)
                    self.r33[i][j][k] = self.r_zeta2(i, j, k)
                    self.r12[i][j][k] = self.r_xi_eta(i, j, k)
                    self.r23[i][j][k] = self.r_eta_zeta(i, j, k)
                    self.r31[i][j][k] = self.r_zeta_xi(i, j, k)
                    self.a1[i][j][k] = self.alpha1(i, j, k)
                    self.a2[i][j][k] = self.alpha2(i, j, k)
                    self.a3[i][j][k] = self.alpha3(i, j, k)
                    self.b12[i][j][k] = self.beta12(i, j, k)
                    self.b23[i][j][k] = self.beta23(i, j, k)
                    self.b31[i][j][k] = self.beta31(i, j, k)

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

    def smooth(self):
        """
        Smooth the grid using Picard Iteration.
        :return: None.
        """

        '''Temp vars'''
        var_num = (self.i_dim - 1) * (self.j_dim - 1) * (self.k_dim - 1)
        coefficient_matrix = np.zeros((var_num, 19))
        rhs = np.zeros((var_num, 3))

        '''Solve the grid iteratively'''
        iteration_cnt = 0
        residual = sys.float_info.max
        while not math.isclose(residual, 0, abs_tol=1e-5):
            '''Calculate all coefficients'''
            self.calc_all_param()

            '''Build Ar=b'''
            coefficient_matrix.fill(0.0)
            rhs.fill(0.0)
            eqn_idx = 0
            for i in range(1, self.i_dim - 1):
                for j in range(1, self.j_dim - 1):
                    for k in range(1, self.k_dim - 1):
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
                                coefficient_matrix[eqn_idx][t] = ca[t]

                        '''Update Equation counter'''
                        eqn_idx += 1

            '''Solve r'''
            u = dsolve.spsolve(coefficient_matrix, rhs)

            '''Update grid and calculate residual'''
            eqn_idx = 0
            for i in range(1, self.i_dim - 1):
                for j in range(1, self.j_dim - 1):
                    for k in range(1, self.k_dim - 1):
                        residual = max(residual, np.linalg.norm(u[eqn_idx] - self.r[i][j][k], np.inf))
                        self.r[i][j][k] = u[eqn_idx]
                        eqn_idx += 1

            iteration_cnt += 1


class ThomasMiddlecoff3D(EllipticGrid3D):
    def __init__(self, grid):
        super(ThomasMiddlecoff3D, self).__init__(grid)

        self.phi = np.zeros((self.i_dim, self.j_dim, self.k_dim))
        self.psi = np.zeros((self.i_dim, self.j_dim, self.k_dim))
        self.omega = np.zeros((self.i_dim, self.j_dim, self.k_dim))

    def smooth(self):
        pass


class EllipticGridTester(unittest.TestCase):
    def test_2d_laplace(self):
        pass

    def test_2d_tm(self):
        pass

    def test_3d_laplace(self):
        pass

    def test_3d_tm(self):
        pass


if __name__ == '__main__':
    unittest.main()