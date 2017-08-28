import unittest
import numpy as np
from math import sin, cos, pi, exp, sqrt, e
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.opt.ga import RealCodedGA


class Rastrigin(object):
    def __call__(self, x):
        n = len(x)
        A = 10
        return A * n + sum(list(map(lambda u: u ** 2 - A * cos(2 * pi * u), x)))

    def show(self):
        x1 = np.linspace(-5.12, 5.12, 100)
        x2 = np.linspace(-5.12, 5.12, 100)

        X, Y = np.meshgrid(x1, x2, indexing='ij')

        val = np.empty((100, 100))
        for i in range(100):
            for j in range(100):
                val[i][j] = self.__call__((x1[i], x2[j]))

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot_surface(X, Y, val, rstride=1, cstride=1, cmap=plt.cm.hot)
        plt.show()


class Ackley(object):
    def __call__(self, x):
        return -20 * exp(-0.2 * sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))) - exp(0.5 * (cos(2 * pi * x[0]) + cos(2 * pi * x[1]))) + e + 20

    def show(self):
        x1 = np.linspace(-5.12, 5.12, 100)
        x2 = np.linspace(-5.12, 5.12, 100)

        X, Y = np.meshgrid(x1, x2, indexing='ij')

        val = np.empty((100, 100))
        for i in range(100):
            for j in range(100):
                val[i][j] = self.__call__((x1[i], x2[j]))

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot_surface(X, Y, val, rstride=1, cstride=1, cmap=plt.cm.hot)
        plt.show()


class Sphere(object):
    def __call__(self, x):
        return sum(list(map(lambda u: u ** 2, x)))


class Rosenbrock(object):
    def __call__(self, x):
        n = len(x)
        if n < 2:
            raise AssertionError("Insufficient input parameters.")

        ret = 0
        for i in range(n - 1):
            ret += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (x[i] - 1) ** 2

        return ret


class RealCodedGATester(unittest.TestCase):
    @staticmethod
    def test_rastrigin():
        rastrigin = Rastrigin()
        f_rastrigin = rastrigin.__call__
        rg = np.array([[-5.12, 5.12],
                       [-5.12, 5.12],
                       [-5.12, 5.12],
                       [-5.12, 5.12]])
        print('Testing Rastrigin function with {} variables ...'.format(len(rg)))
        rga = RealCodedGA(rg, f_rastrigin, lambda u: -f_rastrigin(u))
        ans = rga.find_optimal(600, 100, 0.05)
        print('Global Minimum: {}, Param: {}\n'.format(ans.value, rga.param_transform(ans.param)))

    @staticmethod
    def test_ackley():
        ackley = Ackley()
        f_ackley = ackley.__call__
        rg = np.array([[-5.12, 5.12],
                       [-5.12, 5.12]])
        print('Testing Ackley function ...')
        rga = RealCodedGA(rg, f_ackley, lambda u: -f_ackley(u))
        ans = rga.find_optimal(300, 60, 0.05)
        print('Global Minimum: {}, Param: {}\n'.format(ans.value, rga.param_transform(ans.param)))

    @staticmethod
    def test_sphere():
        sphere = Sphere()
        f_sphere = sphere.__call__
        rg = np.array([[-5.12, 5.12],
                       [-5.12, 5.12]])
        print('Testing Sphere function with {} variables ...'.format(len(rg)))
        rga = RealCodedGA(rg, f_sphere, lambda u: -f_sphere(u))
        ans = rga.find_optimal(300, 60, 0.05)
        print('Global Minimum: {}, Param: {}\n'.format(ans.value, rga.param_transform(ans.param)))

    @staticmethod
    def test_rosenbrock():
        rosenbrock = Rosenbrock()
        f_rosenbrock = rosenbrock.__call__
        rg = np.array([[-5.12, 5.12],
                       [-5.12, 5.12],
                       [-5.12, 5.12]])
        print('Testing Rosenbrock function with {} variables ...'.format(len(rg)))
        rga = RealCodedGA(rg, f_rosenbrock, lambda u: -f_rosenbrock(u))
        ans = rga.find_optimal(1000, 300, 0.42)
        print('Global Minimum: {}, Param: {}\n'.format(ans.value, rga.param_transform(ans.param)))


if __name__ == '__main__':
    unittest.main()
