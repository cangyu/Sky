import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import dsolve
from src.msh.linear_tfi import Linear_TFI_2D


class Laplace_2D(object):
    def __init__(self, c1, c2, c3, c4, N: int, M: int, zeta=1.0, eta=1.0):
        self.zeta = zeta
        self.eta = eta
        self.N = N
        self.M = M
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4

        '''Initialization'''
        u_list = np.linspace(0, 1.0, M + 1)
        v_list = np.linspace(0, 1.0, N + 1)
        self.r = Linear_TFI_2D(c1, c2, c3, c4).calc_msh(u_list, v_list)

    def calc_msh(self):
        get_diff = lambda a, b: math.pow(a - b, 2)
        residual = 1.0
        while residual > 1e-8:
            cu = self.iterate()
            cnt = 0
            residual = 0.0
            for i in range(1, self.N):
                for j in range(1, self.M):
                    residual += get_diff(cu[0][cnt], self.r[i][j][0])
                    residual += get_diff(cu[1][cnt], self.r[i][j][1])
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
        A = sparse.lil_matrix((unknown_num, unknown_num))
        b = np.zeros(2, unknown_num)

        get_coef = lambda i, j: np.array([-2 * (alpha[i][j] / math.pow(self.zeta, 2) + gamma[i][j] / math.pow(self.eta, 2)),
                                          gamma[i][j] / math.pow(self.eta, 2),
                                          alpha[i][j] / math.pow(self.zeta, 2),
                                          gamma[i][j] / math.pow(self.eta, 2),
                                          alpha[i][j] / math.pow(self.zeta, 2),
                                          -beta[i][j] / 2 / self.eta / self.zeta,
                                          beta[i][j] / 2 / self.zeta / self.beta,
                                          -beta[i][j] / 2 / self.eta / self.zeta,
                                          beta[i][j] / 2 / self.eta / self.zeta])
        get_index_list = lambda k: np.array([k,
                                             k + 1,
                                             k + self.M - 1,
                                             k - 1,
                                             k - self.M + 1,
                                             k + self.M,
                                             k + self.M - 2,
                                             k - self.M,
                                             k - self.M + 2])
        get_pos = lambda i: (int(i / (self.M - 1)) + 1, i % (self.M - 1) + 1)
        is_special = lambda row, col: True if row == 0 or row == self.N or col == 0 or col == self.M else False
        ck = 0
        for i in range(1, self.N):
            for j in range(1, self.M):
                a = get_coef(i, j)
                index = get_index_list(i, j)
                for k in range(0, 9):
                    row, col = get_pos(index[k])
                    if is_special(row, col):
                        b[0][ck] -= a[k] * self.r[row][col][0]
                        b[1][ck] -= a[k] * self.r[row][col][1]
                    else:
                        A[ck][index[k]] = a[k]
            ck += 1

        u = np.zeros(2, unknown_num)
        u[0] = dsolve.spsolve(A, b[0], use_umfpack=True)
        u[1] = dsolve.spsolve(A, b[1], use_umfpack=True)

        return u
