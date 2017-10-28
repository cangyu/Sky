import sys
import numpy as np
from abc import abstractmethod
import pyamg
from scipy import sparse
from scipy.sparse.linalg import dsolve
import math

'''
统一约定：
1. (i,j,k)对应(x,y,z),(u,v.w),(xi, eta, zeta),(x1,x2,x3), (I,J,K)...
2. 不用(N,M),避免行列混淆
'''


class EllipticGrid2D(object):
    def __init__(self, grid):
        """
        二维曲边结构网格
        :param grid: 二维初始网格点数组，下标依次迭代I, J. 每个元素包含(X, Y, Z) / (X, Y)
        """

        if len(grid.shape) != 3:
            raise AssertionError("Invalid input grid!")

        self.I, self.J, self.Dim = grid.shape
        self.r = np.zeros((self.I, self.J, 3))  # 不论输入数组中每个元素是(X, Y)还是(X, Y, Z)，均转换成3维存储
        for i in range(self.I):
            for j in range(self.J):
                for d in range(self.Dim):
                    self.r[i][j][d] = grid[i][j][d]

        self.alpha = np.zeros((self.I, self.J))
        self.beta = np.zeros((self.I, self.J))
        self.gamma = np.zeros((self.I, self.J))

        self.I -= 1  # 网格点在第1维度的最大下标
        self.J -= 1  # 网格点在第2维度的最大下标

        self.delta_xi = 1.0
        self.delta_eta = 1.0
        self.eps = 1e-5  # 残差限

    def set_step(self, xi, eta):
        self.delta_xi = xi
        self.delta_eta = eta

    def set_eps(self, eps):
        self.eps = eps

    def get_grid(self):
        return self.r

    @abstractmethod
    def coefficient_around(self, i, j):
        pass

    @classmethod
    def coordinate_around(cls, i, j):
        """
        (i,j)周围的9个点的坐标
        """

        return np.array([(i, j), (i + 1, j), (i, j + 1), (i - 1, j), (i, j - 1), (i + 1, j + 1), (i - 1, j + 1), (i - 1, j - 1), (i + 1, j - 1)])

    def index_around(self, k):
        """
        第k个点周围9个点的序号 
        """

        return np.array([k, k + 1, k + self.I - 1, k - 1, k - self.I + 1, k + self.I, k + self.I - 2, k - self.I, k - self.I + 2])

    def is_special(self, i, j):
        """
        判断坐标是否在边界上
        """

        return True if i in (0, self.I) or j in (0, self.J) else False

    def pder(self, i, j, comp, dir1, dir2=None):
        """
        中心差分偏导数
        :param i: 网格点在第1维度的下标
        :param j: 网格点在第2维度的下标
        :param comp: 分量
        :param dir1: 偏导数方向
        :param dir2: 偏导数方向
        :return: 在(i,j)处的偏导数
        """

        if comp not in ('x', 'y'):
            raise AssertionError("Invalid component!")
        if dir1 not in ('xi', 'eta'):
            raise AssertionError("Invalid first direction!")
        if (dir2 is not None) and (dir2 not in ('xi', 'eta')):
            raise AssertionError("Invalid second direction!")

        comp = ord(comp) - ord('x')

        if dir2 is None:
            if dir1 == 'xi':
                return (self.r[i + 1][j][comp] - self.r[i - 1][j][comp]) / (2.0 * self.delta_xi)
            else:
                return (self.r[i][j + 1][comp] - self.r[i][j - 1][comp]) / (2.0 * self.delta_eta)
        else:
            if dir1 == dir2:
                if dir1 == 'xi':
                    return (self.r[i + 1][j][comp] - 2.0 * self.r[i][j][comp] + self.r[i - 1][j][comp]) / (self.delta_xi ** 2)
                else:
                    return (self.r[i][j + 1][comp] - 2.0 * self.r[i][j][comp] + self.r[i][j - 1][comp]) / (self.delta_eta ** 2)
            else:
                return (self.r[i + 1][j + 1][comp] + self.r[i - 1][j - 1][comp] - self.r[i + 1][j - 1][comp] - self.r[i - 1][j + 1][comp]) / (4.0 * self.delta_xi * self.delta_eta)

    def pvec(self, i, j, dir1, dir2=None):
        """
        中心差分偏导矢
        :param i: 网格点在第1维度的下标
        :param j: 网格点在第2维度的下标
        :param dir1: 偏导数方向
        :param dir2: 偏导数方向
        :return: 在(i,j)处的偏导数
        """

        if dir1 not in ('xi', 'eta'):
            raise AssertionError("Invalid first direction!")
        if (dir2 is not None) and (dir2 not in ('xi', 'eta')):
            raise AssertionError("Invalid second direction!")

        if dir2 is None:
            if dir1 == 'xi':
                return (self.r[i + 1][j] - self.r[i - 1][j]) / (2.0 * self.delta_xi)
            else:
                return (self.r[i][j + 1] - self.r[i][j - 1]) / (2.0 * self.delta_eta)
        else:
            if dir1 == dir2:
                if dir1 == 'xi':
                    return (self.r[i + 1][j] - 2.0 * self.r[i][j] + self.r[i - 1][j]) / (self.delta_xi ** 2)
                else:
                    return (self.r[i][j + 1] - 2.0 * self.r[i][j] + self.r[i][j - 1]) / (self.delta_eta ** 2)
            else:
                return (self.r[i + 1][j + 1] - self.r[i + 1][j - 1] - self.r[i - 1][j + 1] + self.r[i - 1][j - 1]) / (4.0 * self.delta_xi * self.delta_eta)

    def calc_alpha_beta_gamma(self):
        """
        计算所有网格点上的alpha，beta, gamma
        """

        for i in range(1, self.I):
            for j in range(1, self.J):
                xx = self.pder(i, j, 'x', 'xi')
                xe = self.pder(i, j, 'x', 'eta')
                yx = self.pder(i, j, 'y', 'xi')
                ye = self.pder(i, j, 'y', 'eta')
                self.alpha[i][j] = xe ** 2 + ye ** 2
                self.beta[i][j] = xx * xe + yx * ye
                self.gamma[i][j] = xx ** 2 + yx ** 2

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

    def calc_grid(self, choice='Picard', w=1.0, use_amg=False):
        """
        SOR迭代次数较多，求解较慢
        Picard迭代基于9点格式，通过求解稀疏矩阵的方式来求解网格，求解较快
        :return: None
        """

        if choice == 'SOR':
            if w is None:
                raise ValueError('Should specify relax factor!')
            else:
                self.sor_iteration(w)
        elif choice == 'Picard':
            if use_amg is None:
                raise ValueError('Should specify whether using a AMG solver！')
            else:
                self.picard_iteration(use_amg)
        else:
            raise ValueError('Invalid solution choice!')

    @abstractmethod
    def sor_iteration(self, w):
        pass

    @abstractmethod
    def picard_iteration(self, use_amg):
        pass

    @classmethod
    def calc_diff(cls, grid1, grid2):
        """
        计算两套网格间差值 
        """

        if grid1.shape != grid2.shape:
            raise AssertionError("Grid shape do not match!")

        tmp = grid1 - grid2
        ans = 0.0
        for i in range(len(tmp)):
            for j in range(len(tmp[i])):
                ans = max(ans, np.linalg.norm(tmp[i][j], np.inf))

        return ans


class Laplace2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        对网格做Laplace光顺
        :param grid: 初始网格
        """

        super(Laplace2D, self).__init__(grid)

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
        while residual > self.eps:
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
    def __init__(self, grid):
        self.r = np.copy(grid)
        self.d_xi = 1.0
        self.d_eta = 1.0
        self.d_zeta = 1.0

    def r_xi(self, i, j, k):
        return 0.5 * (self.r[i + 1][j][k] - self.r[i - 1][j][k]) / self.d_xi

    def r_eta(self, i, j, k):
        return 0.5 * (self.r[i][j + 1][k] - self.r[i][j - 1][k]) / self.d_eta

    def r_zeta(self, i, j, k):
        return 0.5 * (self.r[i][j][k + 1] - self.r[i][j][k - 1]) / self.d_zeta

    def r_xi2(self, i, j, k):
        return (self.r[i - 1][j][k] - 2 * self.r[i][j][k] + self.r[i + 1][j][k]) / math.pow(self.d_xi, 2)

    def r_eta2(self, i, j, k):
        return (self.r[i][j - 1][k] - 2 * self.r[i][j][k] + self.r[i][j + 1][k]) / math.pow(self.d_eta, 2)

    def r_zeta2(self, i, j, k):
        return (self.r[i][j][k - 1] - 2 * self.r[i][j][k] + self.r[i][j][k + 1]) / math.pow(self.d_zeta, 2)

    def r_xi_eta(self, i, j, k):
        return 0.25 * (self.r[i + 1][j + 1][k] - self.r[i + 1][j - 1][k] - self.r[i - 1][j + 1][k] + self.r[i - 1][j - 1][k]) / (self.d_xi * self.d_eta)

    def r_eta_zeta(self, i, j, k):
        return 0.25 * (self.r[i][j + 1][k + 1] - self.r[i][j - 1][k + 1] - self.r[i][j + 1][k - 1] + self.r[i][j - 1][k - 1]) / (self.d_eta * self.d_zeta)

    def r_zeta_xi(self, i, j, k):
        return 0.25 * (self.r[i + 1][j][k + 1] - self.r[i - 1][j][k + 1] - self.r[i + 1][j][k - 1] + self.r[i - 1][j][k - 1]) / (self.d_zeta * self.d_xi)

    def alpha1(self, i, j, k):
        pass

    def alpha2(self, i, j, k):
        pass

    def alpha3(self, i, j, k):
        pass

    def beta12(self, i, j, k):
        pass

    def beta23(self, i, j, k):
        pass

    def beta31(self, i, j, k):
        pass

    def jacobi(self, i, j, k):
        t = np.matrix()
