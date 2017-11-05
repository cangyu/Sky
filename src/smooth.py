import unittest
import sys
import math
import numpy as np
from abc import abstractmethod
from scipy.linalg import norm
from scipy import sparse
from scipy.sparse.linalg import dsolve

"""
Implementation of the grid smoothing tools using elliptic PDE operator.

Note:
1. (i,j,k) is corresponding to (x,y,z),(u,v.w),(xi, eta, zeta),(x1,x2,x3),(I,J,K)...
2. (N,M) notation is not suggested to avoid confuse on column-and-row.
3. By default, delta_xi = delta_eta = delta_zeta = 1
4. All the partial derivatives in elliptic PDE are calculated with the central scheme.
"""


def vector_square(u):
    return sum(map(lambda x: x ** 2, u))


class EllipticGrid2D(object):
    di = [0, 1, 0, -1, 0, 1, -1, -1, 1]
    dj = [0, 0, 1, 0, -1, 1, 1, -1, -1]

    def __init__(self, grid):
        """
        二维曲边结构网格
        :param grid: 二维初始网格点数组，下标依次迭代I, J. 每个元素包含(X, Y, Z) / (X, Y)
        """

        '''Pre-check'''
        if len(grid.shape) != 3 or grid.shape[-1] > 3:
            raise AssertionError("Invalid input grid!")

        '''Shape constants'''
        ii, jj, dim = grid.shape

        '''Grid'''
        self.r = np.copy(grid)

        '''Partial derivatives'''
        self.r1 = np.empty_like(grid)
        self.r2 = np.empty_like(grid)
        self.r11 = np.empty_like(grid)
        self.r22 = np.empty_like(grid)
        self.r12 = np.empty_like(grid)

        '''Source term'''
        self.pq = np.empty_like(grid)

        '''Coefficients'''
        self.alpha = np.zeros((ii, jj))
        self.beta = np.zeros((ii, jj))
        self.gamma = np.zeros((ii, jj))
        self.jacobi = np.zeros((ii, jj))

    @property
    def i_dim(self):
        return self.r.shape[0]

    @property
    def j_dim(self):
        return self.r.shape[1]

    def is_special(self, i, j):
        return i == 0 or i == self.i_dim - 1 or j == 0 or j == self.j_dim - 1

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

    def alpha(self, i, j, k):
        return vector_square(self.r2[i][j][k])

    def beta(self, i, j, k):
        return np.dot(self.r2[i][j][k], self.r1[i][j][k])

    def gamma(self, i, j, k):
        return vector_square(self.r1[i][j][k])

    def calc_coefficient_matrix(self):
        """
        Calculate the coefficient matrix for the 9 point stencil
        """

        var_cnt = (self.I - 1) * (self.J - 1)
        b = np.zeros((2, var_cnt))
        rows = []
        cols = []
        val = []
        ck = 0  # equation counter

        for j in range(1, self.J):
            for i in range(1, self.I):
                a = self.coefficient_around(i, j)
                coord = self.coordinate_around(i, j)
                index = self.index_around(ck)
                for k in range(9):
                    ii, jj = coord[k]
                    if self.is_special(ii, jj):
                        b[0][ck] -= a[k] * self.r[ii][jj][0]
                        b[1][ck] -= a[k] * self.r[ii][jj][1]
                    else:
                        rows.append(ck)
                        cols.append(index[k])
                        val.append(a[k])

                ck += 1

        '''Construct sparse coefficient matrix'''
        scm = sparse.coo_matrix((val, (rows, cols)), shape=(var_cnt, var_cnt), dtype=float)
        scm_csr = scm.tocsr()
        return scm_csr, b

    @abstractmethod
    def smooth(self):
        pass


class Laplace2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        对网格做Laplace光顺
        :param grid: 初始网格
        """

        super(Laplace2D, self).__init__(grid)

    def smooth(self):


    def sor_iteration(self, w):
        """
        超松弛迭代(SOR)计算内部网格点
        :param w: 松弛因子
        :return: None
        """

        iteration_cnt = 0
        residual = sys.float_info.max
        while residual > self.eps:
            self.calc_alpha_beta_gamma()
            new_grid = np.copy(self.r)

            for i in range(1, self.I):
                for j in range(1, self.J):
                    t1 = self.alpha[i][j] / (self.delta_xi ** 2) * (self.r[i + 1][j] + self.r[i - 1][j])
                    t2 = 2.0 * self.beta[i][j] * np.array([self.pder(i, j, 'x', 'xi', 'eta'), self.pder(i, j, 'y', 'xi', 'eta'), 0])
                    t3 = self.gamma[i][j] / (self.delta_eta ** 2) * (self.r[i][j + 1] + self.r[i][j - 1])
                    t0 = 2.0 * (self.alpha[i][j] / (self.delta_xi ** 2) + self.gamma[i][j] / (self.delta_eta ** 2))
                    t = (t1 - t2 + t3) / t0
                    new_grid[i][j] = w * t + (1.0 - w) * self.r[i][j]

            iteration_cnt += 1
            residual = EllipticGrid2D.calc_diff(self.r, new_grid)
            print("{}: {}".format(iteration_cnt, residual))
            self.r = np.copy(new_grid)

    def coefficient_around(self, i, j):
        ax2 = self.alpha[i][j] / self.delta_xi ** 2
        ge2 = self.gamma[i][j] / self.delta_eta ** 2
        bxe2 = self.beta[i][j] / (2.0 * self.delta_xi * self.delta_eta)

        return np.array([-2 * (ax2 + ge2), ax2, ge2, ax2, ge2, -bxe2, bxe2, -bxe2, bxe2], float)

    def picard_iteration(self, use_amg):
        iteration_cnt = 0
        residual = sys.float_info.max
        while residual > self.eps:
            u = np.zeros((2, (self.I - 1) * (self.J - 1)))
            old_grid = np.copy(self.r)

            self.calc_alpha_beta_gamma()
            A, b = self.calc_coefficient_matrix()

            u[0] = dsolve.spsolve(A, b[0], use_umfpack=True)
            u[1] = dsolve.spsolve(A, b[1], use_umfpack=True)

            cnt = 0
            for j in range(1, self.J):
                for i in range(1, self.I):
                    self.r[i][j][0] = u[0][cnt]
                    self.r[i][j][1] = u[1][cnt]
                    cnt += 1

            iteration_cnt += 1
            residual = EllipticGrid2D.calc_diff(self.r, old_grid)
            print("{}: {}".format(iteration_cnt, residual))


class ThomasMiddlecoff2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        基于Thomas-Middlecoff方法对网格进行光顺
        """

        super(ThomasMiddlecoff2D, self).__init__(grid)
        self.phi = np.zeros((self.I + 1, self.J + 1))
        self.psi = np.zeros((self.I + 1, self.J + 1))

    def calc_phi_psi(self):
        """
        计算所有网格点上的phi, psi
        """

        '''Boundary'''
        for j in (0, self.J):
            for i in range(1, self.I):
                xx = self.pder(i, j, 'x', 'xi')
                yx = self.pder(i, j, 'y', 'xi')
                xxx = self.pder(i, j, 'x', 'xi', 'xi')
                yxx = self.pder(i, j, 'y', 'xi', 'xi')
                self.phi[i][j] = -(xx * xxx + yx * yxx) / (xx ** 2 + yx ** 2)

        for i in (0, self.I):
            for j in range(1, self.J):
                xe = self.pder(i, j, 'x', 'eta')
                ye = self.pder(i, j, 'y', 'eta')
                xee = self.pder(i, j, 'x', 'eta', 'eta')
                yee = self.pder(i, j, 'y', 'eta', 'eta')
                self.psi[i][j] = -(xe * xee + ye * yee) / (xe ** 2 + ye ** 2)

        '''Linear Interpolate'''
        for i in range(1, self.I):
            dist = np.linspace(self.phi[i][0], self.phi[i][self.J], self.J + 1)
            for j in range(1, self.J):
                self.phi[i][j] = dist[j]

        for j in range(1, self.J):
            dist = np.linspace(self.psi[0][j], self.psi[self.I][j], self.I + 1)
            for i in range(1, self.I):
                self.psi[i][j] = dist[i]

    def sor_iteration(self, w):
        """
        超松弛迭代计算内部网格点
        :param w: 松弛因子
        :return: None
        """

        iteration_cnt = 0
        residual = sys.float_info.max
        while residual > self.eps:
            self.calc_alpha_beta_gamma()
            self.calc_phi_psi()
            new_grid = np.copy(self.r)

            for i in range(1, self.I):
                for j in range(1, self.J):
                    t1 = self.alpha[i][j] * (1.0 / self.delta_xi ** 2 + self.phi[i][j] / (2 * self.delta_xi)) * self.r[i + 1][j]
                    t2 = self.alpha[i][j] * (1.0 / self.delta_xi ** 2 - self.phi[i][j] / (2 * self.delta_xi)) * self.r[i - 1][j]
                    t3 = 2 * self.beta[i][j] * self.pvec(i, j, 'xi', 'eta')
                    t4 = self.gamma[i][j] * (1.0 / self.delta_eta ** 2 + self.psi[i][j] / (2 * self.delta_eta)) * self.r[i][j + 1]
                    t5 = self.gamma[i][j] * (1.0 / self.delta_eta ** 2 - self.psi[i][j] / (2 * self.delta_eta)) * self.r[i][j - 1]
                    t0 = 2 * (self.alpha[i][j] / self.delta_xi ** 2 + self.gamma[i][j] / self.delta_eta ** 2)
                    t = (t1 + t2 - t3 + t4 + t5) / t0
                    new_grid[i][j] = w * t + (1 - w) * self.r[i][j]

            iteration_cnt += 1
            residual = EllipticGrid2D.calc_diff(self.r, new_grid)
            print("{}: {}".format(iteration_cnt, residual))
            self.r = np.copy(new_grid)

    def coefficient_around(self, i, j):
        ax2 = self.alpha[i][j] / self.delta_xi ** 2
        ge2 = self.gamma[i][j] / self.delta_eta ** 2
        bxe2 = self.beta[i][j] / (2.0 * self.delta_xi * self.delta_eta)
        x21 = 1.0 / self.delta_xi ** 2
        e21 = 1.0 / self.delta_eta ** 2
        phx2 = self.phi[i][j] / (2.0 * self.delta_xi)
        pse2 = self.psi[i][j] / (2.0 * self.delta_eta)

        return np.array([-2.0 * (ax2 + ge2),
                         self.alpha[i][j] * (x21 + phx2),
                         self.gamma[i][j] * (e21 + pse2),
                         self.alpha[i][j] * (x21 - phx2),
                         self.gamma[i][j] * (e21 - pse2),
                         -bxe2, bxe2, -bxe2, bxe2])

    def picard_iteration(self, use_amg):
        """
        9点格式矩阵迭代计算内部网格点
        """

        iteration_cnt = 0
        residual = sys.float_info.max
        while not math.isclose(residual, 0, abs_tol=1e-5):
            u = np.zeros((2, (self.I - 1) * (self.J - 1)))
            old_grid = np.copy(self.r)

            self.calc_alpha_beta_gamma()
            self.calc_phi_psi()
            A, b = self.calc_coefficient_matrix()

            if use_amg:
                ml = pyamg.ruge_stuben_solver(A)
                u[0] = ml.solve(b[0], tol=1e-10)
                u[1] = ml.solve(b[1], tol=1e-10)
            else:
                u[0] = dsolve.spsolve(A, b[0], use_umfpack=True)
                u[1] = dsolve.spsolve(A, b[1], use_umfpack=True)

            cnt = 0
            for j in range(1, self.J):
                for i in range(1, self.I):
                    self.r[i][j][0] = u[0][cnt]
                    self.r[i][j][1] = u[1][cnt]
                    cnt += 1

            iteration_cnt += 1
            residual = EllipticGrid2D.calc_diff(self.r, old_grid)
            print("{}: {}".format(iteration_cnt, residual))


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
                        self.j[i][j][k] = self.jacobi(i, j, k)

            '''Build Ar=b'''
            coefficient_matrix.fill(0.0)
            rhs.fill(0.0)
            eqn_idx = 0
            for i in range(1, self.i_dim - 1):
                for j in range(1, self.j_dim - 1):
                    for k in range(1, self.k_dim - 1):
                        '''Calculate stencil coefficients'''
                        ca = np.empty(19, float)
                        ca[0] = -2 * (self.a1[i][j][k] + self.a2[i][j][k] + self.a3[i][j][k])
                        ca[1] = ca[2] = self.a1[i][j][k]
                        ca[3] = ca[4] = self.a2[i][j][k]
                        ca[5] = ca[6] = self.a3[i][j][k]
                        ca[7] = ca[8] = 0.5 * self.b12[i][j][k]
                        ca[9] = ca[10] = -ca[8]
                        ca[11] = ca[12] = 0.5 * self.b23[i][j][k]
                        ca[13] = ca[14] = -ca[12]
                        ca[15] = ca[16] = 0.5 * self.b31[i][j][k]
                        ca[17] = ca[18] = -ca[16]

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
