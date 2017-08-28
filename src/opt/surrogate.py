import numpy as np
from numpy.linalg import inv
import math
from scipy.linalg import norm
from src.opt.ga import RealCodedGA
from numpy.linalg import det
from abc import abstractmethod


class SurrogateModel(object):
    def __init__(self):
        pass

    @abstractmethod
    def interp(self, x0):
        pass


class Kriging(SurrogateModel):
    def __init__(self, x, z):
        """
        Kriging Surrogate Model.
        :param x: Sample points.
        :param z: Values on sample points.
        """

        '''Pre-Check'''
        n = len(x)
        if len(x) != len(z) or n == 0:
            raise AssertionError("Invalid input.")
        nd = len(x[0])
        f = np.ones(n, float)

        super(Kriging, self).__init__()

        '''Sample points and values'''
        self.x = np.copy(x)
        self.z = np.copy(z)

        '''Distance between sample points'''
        d = np.empty((n, n, nd), float)
        for i in range(n):
            for j in range(i, n):
                for k in range(nd):
                    d[i][j][k] = d[j][i][k] = math.fabs(self.x[i][k] - self.x[j][k])

        def mle(_t):
            _r = np.empty((n, n), float)
            for _ci in range(n):
                for _cj in range(_ci, n):
                    _pw = 0
                    for _ck in range(nd):
                        _pw += _t[_ck] * d[_ci][_cj][_ck] ** 2
                    _r[_ci][_cj] = _r[_cj][_ci] = math.exp(-_pw)

            _r_inv = inv(_r)
            _tmp1 = np.dot(f, _r_inv)
            _expect = np.dot(_tmp1, self.z) / np.dot(_tmp1, f)
            _tmp2 = self.z - _expect * f
            _variance = np.dot(np.dot(_tmp2, _r_inv), _tmp2) / n
            return -(n * math.log(_variance) + math.log(det(_r))) / 2

        def pmle(_t):
            _param = np.copy(_t)
            for _k, _p in enumerate(_param):
                _param[_k] = math.pow(10, _p)
            return mle(_param)

        '''Determine theta_k using GA'''
        rg = np.copy([(-3, 1.3)] * nd)
        rga = RealCodedGA(rg, pmle, pmle)
        self.theta = rga.find_optimal(100 * nd, 20 * nd)
        for _k, _p in enumerate(self.theta):
            self.theta[_k] = math.pow(10, _p)

        '''Correlation under Gauss function'''
        self.R = np.empty((n, n), float)
        for i in range(n):
            for j in range(i, n):
                pdx = 0
                for k in range(nd):
                    pdx += self.theta[k] * d[i][j][k] ** 2
                self.R[i][j] = self.R[j][i] = math.exp(-pdx)

        '''Expect and variance'''
        r_inv = inv(self.R)
        tmp1 = np.dot(f, r_inv)
        e = np.dot(tmp1, self.z) / np.dot(tmp1, f)
        tmp2 = self.z - e * f
        self.var = np.inner(np.dot(tmp2, r_inv), tmp2) / n

        '''Covariance'''
        self.cov = self.var * self.R

        '''Semi-Variance'''
        self.r = np.full((n, n), self.var) - self.cov

    def interp(self, x0):
        """
        Calculate the response value at given point.
        :param x0: Observation point.
        :return: Response value.
        """

        n = len(self.x)
        nd = len(x0)

        r0 = np.empty(n, float)
        for i in range(n):
            tmp = 0.
            for d in range(nd):
                tmp += self.theta[d] * (x0[d] - self.x[i][d]) ** 2
            r0[i] = self.var - self.var * math.exp(-tmp)

        '''Matrix Coefficients'''
        A = np.empty((n + 1, n + 1), float)
        for i in range(n):
            for j in range(n):
                A[i][j] = self.r[i][j]
        for i in range(n):
            A[-1][i] = 1
        for j in range(n):
            A[j][-1] = 1
        A[-1][-1] = 0

        '''RHS'''
        b = np.empty(n + 1, float)
        for i in range(n):
            b[i] = r0[i]
        b[-1] = 1

        '''Solve linear system: 'Ax = b' '''
        x = np.dot(b, inv(A.transpose()))

        '''Assemble'''
        ans = 0
        for i in range(n):
            ans += x[i] * self.z[i]

        return ans

    def show(self):
        pass
