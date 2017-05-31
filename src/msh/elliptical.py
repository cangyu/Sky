import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import dsolve
from src.msh.tfi import Linear_TFI_2D
from src.msh.plot3d import Plot3D, Plot3D_SingleBlock


def square_diff(a, b):
    return math.pow(a - b, 2)


class Laplace_2D(object):
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

        self.zeta = zeta
        self.eta = eta
        self.M = len(pu) - 1
        self.N = len(pv) - 1
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4

        '''Initialize'''
        init_msh = Linear_TFI_2D(c1, c2, c3, c4).calc_msh(pu, pv)
        self.r = np.zeros((self.N + 1, self.M + 1, 2))
        for i in range(0, self.N + 1):
            for j in range(0, self.M + 1):
                self.r[i][j] = init_msh[j][i]

    def calc_msh(self):
        residual = 1.0
        while residual > 1e-12:
            cu = self.iterate()
            cnt = 0
            residual = 0.0
            for i in range(1, self.N):
                for j in range(1, self.M):
                    residual += square_diff(cu[0][cnt], self.r[i][j][0])
                    residual += square_diff(cu[1][cnt], self.r[i][j][1])
                    self.r[i][j][0] = cu[0][cnt]
                    self.r[i][j][1] = cu[1][cnt]
                    cnt += 1

    def iterate(self):
        pder_x_zeta = lambda i, j: 0.5 * (self.r[i + 1][j][0] - self.r[i - 1][j][0]) / self.zeta
        pder_x_eta = lambda i, j: 0.5 * (self.r[i][j + 1][0] - self.r[i][j - 1][0]) / self.eta
        pder_y_zeta = lambda i, j: 0.5 * (self.r[i + 1][j][1] - self.r[i - 1][j][1]) / self.zeta
        pder_y_eta = lambda i, j: 0.5 * (self.r[i][j + 1][1] - self.r[i][j - 1][1]) / self.eta

        '''Pre-compute alpha, beta, gama for each node'''
        alpha = np.zeros((self.N - 1, self.M - 1))
        beta = np.zeros((self.N - 1, self.M - 1))
        gamma = np.zeros((self.N - 1, self.M - 1))

        for i in range(0, self.N - 1):
            ci = i + 1
            for j in range(0, self.M - 1):
                cj = j + 1
                alpha[i][j] = math.pow(pder_x_eta(ci, cj), 2) + math.pow(pder_y_eta(ci, cj), 2)
                beta[i][j] = pder_x_zeta(ci, cj) * pder_x_eta(ci, cj) + pder_y_zeta(ci, cj) * pder_y_eta(ci, cj)
                gamma[i][j] = math.pow(pder_x_zeta(ci, cj), 2) + math.pow(pder_y_zeta(ci, cj), 2)

        '''Coefficients of the equation set'''
        unknown_num = (self.N - 1) * (self.M - 1)
        b = np.zeros((2, unknown_num))

        zt2, et2, zet = math.pow(self.zeta, 2), math.pow(self.eta, 2), self.eta * self.zeta
        coef = lambda i, j: np.array([-2 * (alpha[i][j] / zt2 + gamma[i][j] / et2),
                                      gamma[i][j] / et2,
                                      alpha[i][j] / zt2,
                                      gamma[i][j] / et2,
                                      alpha[i][j] / zt2,
                                      -beta[i][j] / 2 / zet,
                                      beta[i][j] / 2 / zet,
                                      -beta[i][j] / 2 / zet,
                                      beta[i][j] / 2 / zet])

        coord_list = lambda i, j: np.array([(i, j),
                                            (i, j + 1),
                                            (i + 1, j),
                                            (i, j - 1),
                                            (i - 1, j),
                                            (i + 1, j + 1),
                                            (i + 1, j - 1),
                                            (i - 1, j - 1),
                                            (i - 1, j + 1)])

        index_list = lambda k: np.array([k,
                                         k + 1,
                                         k + self.M - 1,
                                         k - 1,
                                         k - self.M + 1,
                                         k + self.M,
                                         k + self.M - 2,
                                         k - self.M,
                                         k - self.M + 2])

        is_special = lambda row, col: True if row == 0 or row == self.N or col == 0 or col == self.M else False

        rows = []
        cols = []
        val = []
        ck = 0
        for i in range(1, self.N):
            for j in range(1, self.M):
                a = coef(i - 1, j - 1)
                coord = coord_list(i, j)
                index = index_list(ck)
                for k in range(0, 9):
                    row, col = coord[k]
                    if is_special(row, col):
                        b[0][ck] -= a[k] * self.r[row][col][0]
                        b[1][ck] -= a[k] * self.r[row][col][1]
                    else:
                        rows.append(ck)
                        cols.append(index[k])
                        val.append(a[k])

                ck += 1

        A = sparse.coo_matrix((val, (rows, cols)), shape=(unknown_num, unknown_num), dtype=float)
        AA = A.tocsr()
        u = np.zeros((2, unknown_num))
        u[0] = dsolve.spsolve(AA, b[0], use_umfpack=True)
        u[1] = dsolve.spsolve(AA, b[1], use_umfpack=True)

        return u

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


class Possion_2D(object):
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

        self.zeta = zeta
        self.eta = eta
        self.M = len(pu) - 1
        self.N = len(pv) - 1
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4

        '''Initialize'''
        init_msh = Linear_TFI_2D(c1, c2, c3, c4).calc_msh(pu, pv)
        self.r = np.zeros((self.N + 1, self.M + 1, 2))
        for i in range(0, self.N + 1):
            for j in range(0, self.M + 1):
                self.r[i][j] = init_msh[j][i]

        '''Solve the grid iteratively'''
        residual = 1.0
        while residual > 1e-12:
            cu = self._iterate()
            cnt = 0
            residual = 0.0
            for i in range(1, self.N):
                for j in range(1, self.M):
                    residual += square_diff(cu[0][cnt], self.r[i][j][0])
                    residual += square_diff(cu[1][cnt], self.r[i][j][1])
                    self.r[i][j][0] = cu[0][cnt]
                    self.r[i][j][1] = cu[1][cnt]
                    cnt += 1

    def _iterate(self):
        """
        迭代计算内部网格点
        :return: 内部网格点坐标
        """

        '''Helpers'''
        pxz = lambda i, j: (self.r[i + 1][j][0] - self.r[i - 1][j][0]) / (2.0 * self.zeta)
        pxe = lambda i, j: (self.r[i][j + 1][0] - self.r[i][j - 1][0]) / (2.0 * self.eta)
        pxzz = lambda i, j: (self.r[i + 1][j][0] - 2 * self.r[i][j][0] + self.r[i - 1][j][0]) / math.pow(self.zeta, 2)
        pxee = lambda i, j: (self.r[i][j + 1][0] - 2 * self.r[i][j][0] + self.r[i][j - 1][0]) / math.pow(self.eta, 2)
        pxze = lambda i, j: (self.r[i + 1][j + 1][0] - self.r[i + 1][j - 1][0] -
                             self.r[i - 1][j + 1][0] + self.r[i - 1][j - 1][0]) / (4 * self.zeta * self.eta)

        pyz = lambda i, j: (self.r[i + 1][j][1] - self.r[i - 1][j][1]) / (2.0 * self.zeta)
        pye = lambda i, j: (self.r[i][j + 1][1] - self.r[i][j - 1][1]) / (2.0 * self.eta)
        pyzz = lambda i, j: (self.r[i + 1][j][1] - 2 * self.r[i][j][1] + self.r[i - 1][j][1]) / math.pow(self.zeta, 2)
        pyee = lambda i, j: (self.r[i][j + 1][1] - 2 * self.r[i][j][1] + self.r[i][j - 1][1]) / math.pow(self.eta, 2)
        pyze = lambda i, j: (self.r[i + 1][j + 1][1] - self.r[i + 1][j - 1][1] -
                             self.r[i - 1][j + 1][1] + self.r[i - 1][j - 1][1]) / (4 * self.zeta * self.eta)

        '''phi, psi'''
        phi = np.zeros((self.N + 1, self.M + 1))
        psi = np.zeros((self.N + 1, self.M + 1))
        cphi = lambda i, j: -(pyz(i, j) * pyzz(i, j) + pxz(i, j) * pxzz(i, j)) / \
                            (math.pow(pxz(i, j), 2) + math.pow(pyz(i, j), 2))
        cpsi = lambda i, j: -(pye(i, j) * pyee(i, j) + pxe(i, j) * pxee(i, j)) / \
                            (math.pow(pxe(i, j), 2) + math.pow(pye(i, j), 2))

        for i in [0, self.N]:
            for j in range(1, self.M):
                phi[i][j] = cphi(i, j)

        for j in [0, self.M]:
            for i in range(1, self.N):
                psi[j][i] = cpsi(j, i)

        for j in range(1, self.M):
            dist = np.linspace(phi[0][j], phi[self.N][j], self.N + 1)
            for i in range(1, self.N):
                phi[i][j] = dist[i]

        for i in range(1, self.N):
            dist = np.linspace(psi[i][0], psi[i][self.M], self.M + 1)
            for j in range(1, self.M):
                psi[i][j] = dist[j]

        '''alpha, beta, gamma'''
        alpha = np.zeros((self.N - 1, self.M - 1))
        beta = np.zeros((self.N - 1, self.M - 1))
        gamma = np.zeros((self.N - 1, self.M - 1))

        for i in range(0, self.N - 1):
            ci = i + 1
            for j in range(0, self.M - 1):
                cj = j + 1
                alpha[i][j] = math.pow(pxe(ci, cj), 2) + math.pow(pye(ci, cj), 2)
                beta[i][j] = pxe(ci, cj) * pxz(ci, cj) + pye(ci, cj) * pyz(ci, cj)
                gamma[i][j] = math.pow(pxz(ci, cj), 2) + math.pow(pyz(ci, cj), 2)

        '''Coefficients of the equation set'''
        unknown_num = (self.N - 1) * (self.M - 1)
        b = np.zeros((2, unknown_num))

        zt2, et2, zet = math.pow(self.zeta, 2), math.pow(self.eta, 2), self.eta * self.zeta
        coef = lambda i, j: np.array([-2 * (alpha[i][j] / zt2 + gamma[i][j] / et2),
                                      gamma[i][j] * (1 / et2 + psi[i][j] / 2 / self.eta),
                                      alpha[i][j] * (1 / zt2 + phi[i][j] / 2 / self.zeta),
                                      gamma[i][j] / (1 / et2 - psi[i][j] / 2 / self.eta),
                                      alpha[i][j] / (1 / zt2 - phi[i][j] / 2 / self.zeta),
                                      -beta[i][j] / 2 / zet,
                                      beta[i][j] / 2 / zet,
                                      -beta[i][j] / 2 / zet,
                                      beta[i][j] / 2 / zet])

        coord_list = lambda i, j: np.array([(i, j),
                                            (i, j + 1),
                                            (i + 1, j),
                                            (i, j - 1),
                                            (i - 1, j),
                                            (i + 1, j + 1),
                                            (i + 1, j - 1),
                                            (i - 1, j - 1),
                                            (i - 1, j + 1)])

        index_list = lambda k: np.array([k,
                                         k + 1,
                                         k + self.M - 1,
                                         k - 1,
                                         k - self.M + 1,
                                         k + self.M,
                                         k + self.M - 2,
                                         k - self.M,
                                         k - self.M + 2])

        is_special = lambda row, col: True if row == 0 or row == self.N or col == 0 or col == self.M else False

        rows = []
        cols = []
        val = []
        ck = 0
        for i in range(1, self.N):
            for j in range(1, self.M):
                a = coef(i - 1, j - 1)
                coord = coord_list(i, j)
                index = index_list(ck)
                for k in range(0, 9):
                    row, col = coord[k]
                    if is_special(row, col):
                        b[0][ck] -= a[k] * self.r[row][col][0]
                        b[1][ck] -= a[k] * self.r[row][col][1]
                    else:
                        rows.append(ck)
                        cols.append(index[k])
                        val.append(a[k])

                ck += 1

        A = sparse.coo_matrix((val, (rows, cols)), shape=(unknown_num, unknown_num), dtype=float)
        AA = A.tocsr()
        u = np.zeros((2, unknown_num))
        u[0] = dsolve.spsolve(AA, b[0], use_umfpack=True)
        u[1] = dsolve.spsolve(AA, b[1], use_umfpack=True)

        return u
