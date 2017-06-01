import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import dsolve
from src.msh.tfi import Linear_TFI_2D
from src.msh.plot3d import Plot3D, Plot3D_SingleBlock
from operator import mul
from functools import reduce


def square_diff(a: float, b: float):
    return math.pow(a - b, 2)


def coord_list(i: int, j: int):
    """
    (i,j)周围的9个点的坐标
    """

    return np.array([(i, j), (i, j + 1), (i + 1, j), (i, j - 1), (i - 1, j),
                     (i + 1, j + 1), (i + 1, j - 1), (i - 1, j - 1), (i - 1, j + 1)])


def index_list(k: int, M: int):
    """
    第k个点周围9个点的序号 
    """

    return np.array([k, k + 1, k + M - 1, k - 1, k - M + 1, k + M, k + M - 2, k - M, k - M + 2])


def is_special(row, col, N, M):
    return True if row in (0, N) or col in (0, M) else False


def pxz(r, i, j, gp):
    return (r[i + 1][j][0] - r[i - 1][j][0]) / (2 * gp)


def pxe(r, i, j, gp):
    return (r[i][j + 1][0] - r[i][j - 1][0]) / (2 * gp)


def pyz(r, i, j, gp):
    return (r[i + 1][j][1] - r[i - 1][j][1]) / (2 * gp)


def pye(r, i, j, gp):
    return (r[i][j + 1][1] - r[i][j - 1][1]) / (2 * gp)


def pxzz(r, i, j, gp):
    return (r[i + 1][j][0] - 2 * r[i][j][0] + r[i - 1][j][0]) / gp ** 2


def pxee(r, i, j, gp):
    return (r[i][j + 1][0] - 2 * r[i][j][0] + r[i][j - 1][0]) / gp ** 2


def pxze(r, i, j, gp):
    return (r[i + 1][j + 1][0] - r[i + 1][j - 1][0] - r[i - 1][j + 1][0] + r[i - 1][j - 1][0]) / (4 * reduce(mul, gp))


def pyzz(r, i, j, gp):
    return (r[i + 1][j][1] - 2 * r[i][j][1] + r[i - 1][j][1]) / gp ** 2


def pyee(r, i, j, gp):
    return (r[i][j + 1][1] - 2 * r[i][j][1] + r[i][j - 1][1]) / gp ** 2


def pyze(r, i, j, gp):
    return (r[i + 1][j + 1][1] - r[i + 1][j - 1][1] - r[i - 1][j + 1][1] + r[i - 1][j - 1][1]) / (4 * reduce(mul, gp))


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
        self.r = np.zeros((self.N + 1, self.M + 1, 2))

        '''Initialize'''
        init_msh = Linear_TFI_2D(c1, c2, c3, c4).calc_msh(pu, pv)
        for i in range(0, self.N + 1):
            for j in range(0, self.M + 1):
                self.r[i][j] = init_msh[j][i]

        '''Solve iteratively'''
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

        unknown_num = (self.N - 1) * (self.M - 1)
        b = np.zeros((2, unknown_num))
        u = np.zeros((2, unknown_num))

        alpha = np.zeros((self.N + 1, self.M + 1))
        beta = np.zeros((self.N + 1, self.M + 1))
        gamma = np.zeros((self.N + 1, self.M + 1))

        zt2 = math.pow(self.zeta, 2)
        et2 = math.pow(self.eta, 2)
        zet = self.eta * self.zeta

        coef = lambda i, j: np.array([-2 * (alpha[i][j] / zt2 + gamma[i][j] / et2),
                                      gamma[i][j] / et2,
                                      alpha[i][j] / zt2,
                                      gamma[i][j] / et2,
                                      alpha[i][j] / zt2,
                                      -beta[i][j] / 2 / zet,
                                      beta[i][j] / 2 / zet,
                                      -beta[i][j] / 2 / zet,
                                      beta[i][j] / 2 / zet])

        '''alpha, beta, gama'''
        for i in range(1, self.N):
            for j in range(1, self.M):
                vpxe = pxe(self.r, i, j, self.eta)
                vpxz = pxz(self.r, i, j, self.zeta)
                vpye = pye(self.r, i, j, self.eta)
                vpyz = pyz(self.r, i, j, self.zeta)

                alpha[i][j] = vpxe * vpxe + vpye * vpye
                beta[i][j] = vpxe * vpxz + vpye * vpyz
                gamma[i][j] = vpxz * vpxz + vpyz * vpyz

        '''Coefficient Matrix'''
        rows = []
        cols = []
        val = []
        ck = 0
        for i in range(1, self.N):
            for j in range(1, self.M):
                a = coef(i, j)
                coord = coord_list(i, j)
                index = index_list(ck, self.M)
                for k in range(0, 9):
                    row, col = coord[k]
                    if is_special(row, col, self.N, self.M):
                        b[0][ck] -= a[k] * self.r[row][col][0]
                        b[1][ck] -= a[k] * self.r[row][col][1]
                    else:
                        rows.append(ck)
                        cols.append(index[k])
                        val.append(a[k])

                ck += 1

        A = sparse.coo_matrix((val, (rows, cols)), shape=(unknown_num, unknown_num), dtype=float)
        AA = A.tocsr()

        '''Solve'''
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

        zt2 = math.pow(self.zeta, 2)
        et2 = math.pow(self.eta, 2)
        zet = self.eta * self.zeta

        alpha = np.zeros((self.N + 1, self.M + 1))
        beta = np.zeros((self.N + 1, self.M + 1))
        gamma = np.zeros((self.N + 1, self.M + 1))

        phi = np.zeros((self.N + 1, self.M + 1))
        psi = np.zeros((self.N + 1, self.M + 1))

        unknown_num = (self.N - 1) * (self.M - 1)
        b = np.zeros((2, unknown_num))
        u = np.zeros((2, unknown_num))

        coef = lambda i, j: np.array([-2 * (alpha[i][j] / zt2 + gamma[i][j] / et2),
                                      gamma[i][j] * (1 / et2 + psi[i][j] / 2 / self.eta),
                                      alpha[i][j] * (1 / zt2 + phi[i][j] / 2 / self.zeta),
                                      gamma[i][j] / (1 / et2 - psi[i][j] / 2 / self.eta),
                                      alpha[i][j] / (1 / zt2 - phi[i][j] / 2 / self.zeta),
                                      -beta[i][j] / 2 / zet,
                                      beta[i][j] / 2 / zet,
                                      -beta[i][j] / 2 / zet,
                                      beta[i][j] / 2 / zet])

        '''alpha, beta, gamma'''
        for i in range(1, self.N):
            for j in range(1, self.M):
                vpxe = pxe(self.r, i, j, self.eta)
                vpxz = pxz(self.r, i, j, self.zeta)
                vpye = pye(self.r, i, j, self.eta)
                vpyz = pyz(self.r, i, j, self.zeta)

                alpha[i][j] = vpxe * vpxe + vpye * vpye
                beta[i][j] = vpxe * vpxz + vpye * vpyz
                gamma[i][j] = vpxz * vpxz + vpyz * vpyz

        '''phi, psi'''
        for i in [0, self.N]:
            for j in range(1, self.M):
                vpxz = pxz(self.r, i, j, self.zeta)
                vpxzz = pxzz(self.r, i, j, self.zeta)
                vpyz = pyz(self.r, i, j, self.zeta)
                vpyzz = pyzz(self.r, i, j, self.zeta)
                phi[i][j] = -(vpxz * vpxzz + vpyz * vpyzz) / (vpxz ** 2 + vpyz ** 2)

        for j in [0, self.M]:
            for i in range(1, self.N):
                vpxe = pxe(self.r, i, j, self.eta)
                vpxee = pxee(self.r, i, j, self.eta)
                vpye = pye(self.r, i, j, self.eta)
                vpyee = pyee(self.r, i, j, self.eta)
                psi[i][j] = -(vpxe * vpxee + vpye * vpyee) / (vpxe ** 2 + vpye ** 2)

        for j in range(1, self.M):
            dist = np.linspace(phi[0][j], phi[self.N][j], self.N + 1)
            for i in range(1, self.N):
                phi[i][j] = dist[i]

        for i in range(1, self.N):
            dist = np.linspace(psi[i][0], psi[i][self.M], self.M + 1)
            for j in range(1, self.M):
                psi[i][j] = dist[j]

        '''Coefficient Matrix'''
        rows = []
        cols = []
        val = []
        ck = 0
        for i in range(1, self.N):
            for j in range(1, self.M):
                a = coef(i, j)
                coord = coord_list(i, j)
                index = index_list(ck, self.M)
                for k in range(0, 9):
                    row, col = coord[k]
                    if is_special(row, col, self.N, self.M):
                        b[0][ck] -= a[k] * self.r[row][col][0]
                        b[1][ck] -= a[k] * self.r[row][col][1]
                    else:
                        rows.append(ck)
                        cols.append(index[k])
                        val.append(a[k])

                ck += 1

        A = sparse.coo_matrix((val, (rows, cols)), shape=(unknown_num, unknown_num), dtype=float)
        AA = A.tocsr()

        '''Solve'''
        u[0] = dsolve.spsolve(AA, b[0], use_umfpack=True)
        u[1] = dsolve.spsolve(AA, b[1], use_umfpack=True)

        return u
