import numpy as np
from numpy.linalg import inv
import math
from scipy.linalg import norm


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
        n = len(self.x)
        if len(x) != len(z):
            raise AssertionError("Inconsistent input.")

        super(Kriging, self).__init__()

        '''Sample points and values'''
        self.x = np.copy(x)
        self.z = np.copy(z)

        '''Distance between sample points'''
        self.d = np.empty((n, n), float)
        for i in range(n):
            for j in range(i, n):
                self.d[i][j] = self.d[j][i] = norm(self.x[i] - self.x[j], 2)

        '''Correlation under Gauss function'''
        self.theta = 1.0
        self.R = np.empty((n, n), float)
        for i in range(n):
            for j in range(i, n):
                self.R[i][j] = self.R[j][i] = math.exp(-self.theta * math.pow(self.d[i][j], 2))

        '''Expect and variance'''
        f = np.ones(n, float)
        r_inv = inv(self.R)
        tmp1 = np.dot(f, r_inv)
        self.e = np.dot(tmp1, self.z) / np.dot(tmp1, f.transpose())
        tmp2 = self.z - self.e * f
        self.var = np.dot(np.dot(tmp2, r_inv), tmp2.transpose()) / n

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
