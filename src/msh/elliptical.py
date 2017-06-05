import sys
import numpy as np
from src.msh.plot3d import PLOT3D_Block
from abc import abstractmethod

'''
统一约定：
1. (i,j,k)对应(x,y,z),(u,v.w),(xi, eta, mu),(x1,x2,x3), (I,J,K)...
2. 不用(N,M),避免行列混淆
3. 计算域上网格各方向步长均为1
'''


class CurvilinearGrid2D(object):
    def __init__(self, grid):
        """
        二维曲边结构网格
        :param grid: 二维初始网格点数组，下标依次迭代I,J
                     每个元素包含(X, Y, Z) / (X, Y)
        """

        if len(grid.shape) != 3:
            raise AssertionError("Invalid input grid!")

        self.I, self.J, self.Dim = grid.shape
        self.alpha = np.zeros((self.I, self.J))
        self.beta = np.zeros((self.I, self.J))
        self.gamma = np.zeros((self.I, self.J))
        self.r = np.zeros((self.I, self.J, 3))  # 不论输入数组中每个元素是(X, Y)还是(X, Y, Z)，均转换成3维存储
        for i in range(self.I):
            for j in range(self.J):
                for d in range(self.Dim):
                    self.r[i][j][d] = grid[i][j][d]

        self.I -= 1  # 网格点在第1维度的最大下标
        self.J -= 1  # 网格点在第2维度的最大下标

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
                return (self.r[i + 1][j][comp] - self.r[i - 1][j][comp]) / 2
            else:
                return (self.r[i][j + 1][comp] - self.r[i][j - 1][comp]) / 2
        else:
            if dir1 == dir2:
                return self.r[i + 1][j][comp] - 2 * self.r[i][j][comp] + self.r[i - 1][j][comp]
            else:
                return (self.r[i + 1][j + 1][comp] + self.r[i - 1][j - 1][comp] -
                        self.r[i + 1][j - 1][comp] - self.r[i - 1][j + 1][comp]) / 4

    def calc_alpha_beta_gamma(self):
        """
        计算所有网格点上的alpha，beta, gamma
        """

        for j in range(1, self.J):
            for i in range(1, self.I):
                c = np.array([self.pder(i, j, 'x', 'xi'), self.pder(i, j, 'x', 'eta'),
                              self.pder(i, j, 'y', 'xi'), self.pder(i, j, 'y', 'eta')])
                self.alpha[i][j] = c[1] ** 2 + c[3] ** 2
                self.beta[i][j] = c[0] * c[1] + c[2] * c[3]
                self.gamma[i][j] = c[0] ** 2 + c[2] ** 2

    @abstractmethod
    def calc_grid(self):
        pass

    @classmethod
    def calc_diff(cls, grid1, grid2):
        """
        计算两套网格间差值 
        """

        if grid1.shape != grid2.shape:
            raise AssertionError("Grid shape do not match!")

        I, J, Dim = grid1.shape

        ans = 0.0
        for i in range(I):
            for j in range(J):
                for k in range(Dim):
                    ans += np.fabs(grid1[i][j][k] - grid2[i][j][k])

        return ans

    def p3d_blk(self):
        return PLOT3D_Block.build_from_3d(self.r)


class Laplace2D(CurvilinearGrid2D):
    def __init__(self, grid):
        """
        对网格做Laplace光顺
        :param grid: 初始网格
        """
        super(Laplace2D, self).__init__(grid)

    def calc_grid(self, eps=1e-10, w=1.0):
        """
        迭代计算内部网格点
        :param eps: 残差限
        :param w: 松弛因子
        :return: 光顺后的新网格
        """

        residual = sys.float_info.max
        while residual > eps:
            self.calc_alpha_beta_gamma()
            new_grid = np.zeros_like(self.r)
            for i in range(1, self.I):
                for j in range(1, self.J):
                    r_xi_eta = np.array([self.pder(i, j, 'x', 'xi', 'eta'), self.pder(i, j, 'y', 'xi', 'eta')])
                    t = (self.alpha[i][j] * (self.r[i + 1][j] + self.r[i - 1][j]) +
                         self.gamma[i][j] * (self.r[i][j + 1] + self.r[i][j - 1]) -
                         2 * self.beta[i][j] * r_xi_eta) / (2 * (self.alpha[i][j] + self.gamma[i][j]))
                    new_grid[i][j] = w * t + (1 - w) * self.r[i][j]
            residual = CurvilinearGrid2D.calc_diff(self.r, new_grid)
            self.r = np.copy(new_grid)

        return self.r


class ThomasMiddlecoff2D(CurvilinearGrid2D):
    def __init__(self, grid):
        """
        基于Thomas-Middlecoff方法对输入网格进行光顺
        """

        super(ThomasMiddlecoff2D, self).__init__(grid)
        self.phi = np.zeros((self.I + 1, self.J + 1))
        self.psi = np.zeros((self.I + 1, self.J + 1))

    def calc_phi_psi(self):
        for i in (0, self.J):
            for j in range(1, self.I):
                c = np.array([self.pder(i, j, 'x', 'xi'), self.pder(i, j), self.pder(i, j), self.pder(i, j)])
                self.phi[i][j] = -(val_pxz * val_pxzz + val_pyz * val_pyzz) / (val_pxz ** 2 + val_pyz ** 2)

        for j in (0, self.I):
            for i in range(1, self.J):
                val_pxe = self.pxe(i, j)
                val_pxee = self.pxee(i, j)
                val_pye = self.pye(i, j)
                val_pyee = self.pyee(i, j)
                self.psi[i][j] = -(val_pxe * val_pxee + val_pye * val_pyee) / (val_pxe ** 2 + val_pye ** 2)

        for j in range(1, self.M):
            dist = np.linspace(self.phi[0][j], self.phi[self.N][j], self.N + 1)
            for i in range(1, self.N):
                self.phi[i][j] = dist[i]

        for i in range(1, self.N):
            dist = np.linspace(self.psi[i][0], self.psi[i][self.M], self.M + 1)
            for j in range(1, self.M):
                self.psi[i][j] = dist[j]

    def calc_grid(self, eps=1e-10, w=1.0):
        """
        迭代计算内部网格点
        :return: 内部网格点坐标
        :param eps: 残差限
        :param w: 松弛因子
        :return: 光顺后的新网格
        """

        residual = sys.float_info.max
        while residual > eps:
            self.calc_alpha_beta_gamma()
            self.calc_phi_psi()
            new_grid = np.zeros_like(self.r)
            for i in range(1, self.I):
                for j in range(1, self.J):
                    r_xi_eta = np.array([self.pder(i, j, 'x', 'xi', 'eta'), self.pder(i, j, 'y', 'xi', 'eta')])
                    t = (self.alpha[i][j] * (self.r[i + 1][j] + self.r[i - 1][j]) +
                         self.gamma[i][j] * (self.r[i][j + 1] + self.r[i][j - 1]) -
                         2 * self.beta[i][j] * r_xi_eta) / (2 * (self.alpha[i][j] + self.gamma[i][j]))
                    new_grid[i][j] = w * t + (1 - w) * self.r[i][j]
            residual = CurvilinearGrid2D.calc_diff(self.r, new_grid)
            self.r = np.copy(new_grid)

        return self.r
