import numpy as np
from numpy.linalg import inv
import math
from scipy.linalg import norm
from src.opt.ga import RealCodedGA
from numpy.linalg import det


class SurrogateModel(object):
    def __init__(self):
        pass


class Kriging(SurrogateModel):
    def __init__(self, x, z):
        """
        Kriging Surrogate Model.
        :param x: Sample points.
        :param z: Value on sample points.
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
        self.d = np.empty((n, n), float)
        for i in range(n):
            for j in range(i, n):
                self.d[i][j] = self.d[j][i] = norm(self.x[i] - self.x[j], 2)

        def mle(theta):
            R = np.empty((n, n), float)
            for i in range(n):
                for j in range(i, n):
                    pdx = 0
                    for k in range(nd):
                        pdx += theta[k] * self.d[i][j] ** 2
                    R[i][j] = R[j][i] = math.exp(-pdx)

            r_inv = inv(R)
            tmp1 = np.dot(f, r_inv)
            e = np.dot(tmp1, self.z) / np.dot(tmp1, f)
            tmp2 = self.z - e * f
            var = np.dot(np.dot(tmp2, r_inv), tmp2) / n
            return -(n * math.log(var) + math.log(det(R))) / 2

        '''Determine theta_k using GA'''
        rg = np.empty((nd, 2), float)
        for i in range(nd):
            rg[i][0] = -3
            rg[i][1] = 2
        emle = lambda t: mle(np.exp(t))
        rga = RealCodedGA(rg, emle, emle)
        self.theta = rga.find_optimal(100 * nd, 20 * nd, 0.05)

        '''Correlation under Gauss function'''
        self.R = np.empty((n, n), float)
        for i in range(n):
            for j in range(i, n):
                pdx = 0
                for k in range(nd):
                    pdx += self.theta[k] * self.d[i][j] ** 2
                self.R[i][j] = self.R[j][i] = math.exp(-pdx)

        '''Expect and variance'''
        r_inv = inv(self.R)
        tmp1 = np.dot(f, r_inv)
        e = np.dot(tmp1, self.z) / np.dot(tmp1, f)
        tmp2 = self.z - e * f
        self.var = np.dot(np.dot(tmp2, r_inv), tmp2) / n

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
        d0 = np.empty(n, float)
        for i in range(n):
            d0[i] = norm(x0 - self.x[i], 2)

        r0 = np.empty(n, float)
        for i in range(n):
            r0[i] = self.var - self.var * math.exp(-self.theta * math.pow(d0[i], 2))

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
