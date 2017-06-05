import math
import sys
import numpy as np
import pyamg
from scipy import sparse
from scipy.sparse.linalg import dsolve
from src.msh.tfi import Linear_TFI_2D
from abc import ABCMeta, abstractmethod

"""
统一约定：
1. (i,j,k)对应(x,y,z),(u,v.w),(xi, eta, mu),(x1,x2,x3), (I,J,K)...
2. 不用(N,M),避免行列混淆
3. 计算域上网格步长均为1
"""


class CurvilinearGrid2D(object):
    def __init__(self):
        self.r = None
        self.I = int(0)  # 网格点在第1维度的最大下标
        self.J = int(0)  # 网格点在第2维度的最大下标
        self.alpha = None
        self.beta = None
        self.gamma = None

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
                return self.r[i + 1][j][comp] - 2 * self.[i][j][comp] + self.r[i - 1][j][comp]
            else:
                return (self.r[i + 1][j + 1][comp] + self.r[i - 1][j - 1][comp] - self.r[i + 1][j - 1][comp] - self.r[i - 1][j + 1][comp]) / 4

    def is_special(self, i, j):
        """
        判断网格点是否在最外围上
        :param i: 网格点在第1维度的下标
        :param j: 网格点在第2维度的下标
        """

        return True if i in (0, self.I) or j in (0, self.J) else False

    @abstractmethod
    def get_coefficient_list(self, i, j):
        pass

    @classmethod
    def get_coordinate_list(cls, i, j):
        """
        (i,j)周围的9个点的坐标
        """

        return np.array([(i, j),
                         (i + 1, j),
                         (i, j + 1),
                         (i - 1, j),
                         (i, j - 1),
                         (i + 1, j + 1),
                         (i - 1, j + 1),
                         (i - 1, j - 1),
                         (i + 1, j - 1)])

    def get_index_list(self, k):
        """
        第k个点周围9个点的序号 
        """

        return np.array([k, k + 1, k + self.J - 1, k - 1, k - self.J + 1, k + self.J, k + self.J - 2, k - self.J, k - self.J + 2])

    def calc_all_alpha_beta_gamma(self):
        for i in range(1, self.N):
            for j in range(1, self.M):
                val_pxz = self.pxz(i, j)
                val_pxe = self.pxe(i, j)
                val_pyz = self.pyz(i, j)
                val_pye = self.pye(i, j)
                self.alpha[i][j] = val_pxe ** 2 + val_pye ** 2
                self.beta[i][j] = val_pxe * val_pxz + val_pye * val_pyz
                self.gamma[i][j] = val_pxz ** 2 + val_pyz ** 2

    def calc_coefficient_matrix(self):
        b = np.zeros((2, self.unknown_num))
        rows = []
        cols = []
        val = []
        ck = 0
        for i in range(1, self.N):
            for j in range(1, self.M):
                a = self.get_coefficient_list(i, j)
                coord = self.get_coordinate_list(i, j)
                index = self.get_index_list(ck)
                for k in range(0, 9):
                    row, col = coord[k]
                    if self.is_special(row, col):
                        b[0][ck] -= a[k] * self.r[row][col][0]
                        b[1][ck] -= a[k] * self.r[row][col][1]
                    else:
                        rows.append(ck)
                        cols.append(index[k])
                        val.append(a[k])

                ck += 1

        A = sparse.coo_matrix((val, (rows, cols)), shape=(self.unknown_num, self.unknown_num), dtype=float)
        AA = A.tocsr()
        return AA, b


def calc_grid(self, eps=1e-10):
    """
    Solve the grid iteratively
    :param eps: 残差
    :param debug: 调试选项
    :param method: 每步计算方法
    :return: None
    """

    residual = sys.float_info.max
    k = 0
    while residual > eps:
        cu = self.iterate()
        cnt = 0
        k += 1
        residual = 0.0
        for i in range(1, self.N):
            for j in range(1, self.M):
                residual += square_diff(cu[0][cnt], self.r[i][j][0])
                residual += square_diff(cu[1][cnt], self.r[i][j][1])
                self.r[i][j][0] = cu[0][cnt]
                self.r[i][j][1] = cu[1][cnt]
                cnt += 1

        print("{}:{}".format(k, residual))


@abstractmethod
def iterate(self):
    pass


def write_plot3d(self, filename="msh.xyz"):
    K, J, I = 1, self.N + 1, self.M + 1
    pts = np.zeros((K, J, I, 3))

    for k in range(0, K):
        for j in range(0, J):
            for i in range(0, I):
                pts[k][j][i][0] = self.r[j][i][0]
                pts[k][j][i][1] = self.r[j][i][1]

    p3d = Plot3D(I, J, K, pts)
    p3d.output(filename)


def plot3d_blk(self):
    K, J, I = 1, self.N + 1, self.M + 1
    pts = np.zeros((K, J, I, 3))

    for k in range(0, K):
        for j in range(0, J):
            for i in range(0, I):
                pts[k][j][i][0] = self.r[j][i][0]
                pts[k][j][i][1] = self.r[j][i][1]

    return Plot3D_SingleBlock(I, J, K, pts)


class Laplace_2D(CurvilinearGrid2D):
    def __init__(self, c1, c2, c3, c4, pu, pv, zeta=1.0, eta=1.0):
        """
        生成二维Laplace网格
        :param c1: 沿u方向的曲线，靠近u轴
        :param c2: 沿v方向的曲线，靠近v轴
        :param c3: 沿u方向的曲线，远离u轴
        :param c4: 沿v方向的曲线，远离v轴
        :param pu: u方向节点
        :param pv: v方向节点
        :param zeta: u方向网格间距
        :param eta: v方向网格间距
        """

        super(Laplace_2D, self).__init__()

        self.zeta = zeta
        self.eta = eta
        self.M = len(pu) - 1
        self.N = len(pv) - 1
        self.unknown_num = (self.N - 1) * (self.M - 1)
        self.r = np.zeros((self.N + 1, self.M + 1, 2))
        self.alpha = np.zeros((self.N + 1, self.M + 1))
        self.beta = np.zeros((self.N + 1, self.M + 1))
        self.gamma = np.zeros((self.N + 1, self.M + 1))

        '''Initialize'''
        init_msh = Linear_TFI_2D(c1, c2, c3, c4).calc_msh(pu, pv)
        for i in range(0, self.N + 1):
            for j in range(0, self.M + 1):
                self.r[i][j] = init_msh[j][i]

    def get_coefficient_list(self, i, j):
        az2 = self.alpha[i][j] / self.zeta ** 2
        ge2 = self.gamma[i][j] / self.eta ** 2
        bze2 = self.beta[i][j] / (2 * self.zeta * self.eta)

        return np.array([-2 * (az2 + ge2), az2, ge2, az2, ge2, -bze2, bze2, -bze2, bze2])

    def iterate(self):
        """
        迭代计算内部网格点
        """

        u = np.zeros((2, self.unknown_num))

        '''alpha, beta, gamma'''
        self.calc_all_alpha_beta_gamma()

        '''coefficient matrix'''
        A, b = self.calc_coefficient_matrix()

        '''solve'''
        u[0] = dsolve.spsolve(A, b[0], use_umfpack=True)
        u[1] = dsolve.spsolve(A, b[1], use_umfpack=True)

        # ml = pyamg.ruge_stuben_solver(A)
        # u[0] = ml.solve(b[0], tol=1e-10)
        # u[1] = ml.solve(b[1], tol=1e-10)

        return u


class Possion_2D(CurvilinearGrid2D):
    def __init__(self, c1, c2, c3, c4, pu, pv, zeta=1.0, eta=1.0):
        """
        基于Thomas-Middlecoff方法生成结构网格
        :param c1: 沿u方向的曲线，靠近u轴
        :param c2: 沿v方向的曲线，靠近v轴
        :param c3: 沿u方向的曲线，远离u轴
        :param c4: 沿v方向的曲线，远离v轴
        :param pu: u方向节点
        :param pv: v方向节点
        :param zeta: u方向网格间距
        :param eta: v方向网格间距
        """

        super(Possion_2D, self).__init__()

        self.zeta = zeta
        self.eta = eta
        self.M = len(pu) - 1
        self.N = len(pv) - 1
        self.unknown_num = (self.N - 1) * (self.M - 1)
        self.r = np.zeros((self.N + 1, self.M + 1, 2))
        self.alpha = np.zeros((self.N + 1, self.M + 1))
        self.beta = np.zeros((self.N + 1, self.M + 1))
        self.gamma = np.zeros((self.N + 1, self.M + 1))
        self.phi = np.zeros((self.N + 1, self.M + 1))
        self.psi = np.zeros((self.N + 1, self.M + 1))

        '''Laplace Initialize
        init_grid = Laplace_2D(c1, c2, c3, c4, pu, pv, zeta, eta)
        init_grid.calc_grid()
        self.r = np.copy(init_grid.r)
        '''

        '''TFI Initialize'''
        init_msh = Linear_TFI_2D(c1, c2, c3, c4).calc_msh(pu, pv)
        for i in range(0, self.N + 1):
            for j in range(0, self.M + 1):
                self.r[i][j] = init_msh[j][i]

    def get_coefficient_list(self, i, j):
        az2 = self.alpha[i][j] / self.zeta ** 2
        ge2 = self.gamma[i][j] / self.eta ** 2
        bze2 = self.beta[i][j] / (2 * self.zeta * self.eta)
        phz2 = self.phi[i][j] / (2 * self.zeta)
        pse2 = self.psi[i][j] / (2 * self.eta)

        return np.array([-2 * (az2 + ge2), (az2 + phz2), (ge2 + pse2), (az2 - phz2), (ge2 - pse2), -bze2, bze2, -bze2, bze2])

    def iterate(self):
        """
        迭代计算内部网格点
        :return: 内部网格点坐标
        """

        u = np.zeros((2, self.unknown_num))

        '''alpha, beta, gamma'''
        self.calc_all_alpha_beta_gamma()

        '''phi, psi'''
        for i in [0, self.N]:
            for j in range(1, self.M):
                val_pxz = self.pxz(i, j)
                val_pxzz = self.pxzz(i, j)
                val_pyz = self.pyz(i, j)
                val_pyzz = self.pyzz(i, j)
                self.phi[i][j] = -(val_pxz * val_pxzz + val_pyz * val_pyzz) / (val_pxz ** 2 + val_pyz ** 2)

        for j in [0, self.M]:
            for i in range(1, self.N):
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

        '''coefficient matrix'''
        A, b = self.calc_coefficient_matrix()

        '''solve'''
        u[0] = dsolve.spsolve(A, b[0], use_umfpack=True)
        u[1] = dsolve.spsolve(A, b[1], use_umfpack=True)

        # ml = pyamg.ruge_stuben_solver(A)
        # u[0] = ml.solve(b[0], tol=1e-10)
        # u[1] = ml.solve(b[1], tol=1e-10)

        return u
