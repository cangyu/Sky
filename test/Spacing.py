import unittest
import math
import numpy as np
import scipy
from scipy import optimize
import matplotlib.pyplot as plt
from src.grid import single_exponential, double_exponential, hyperbolic_sine, hyperbolic_tangent


class SpacingTestCase(unittest.TestCase):
    def test_single_exponential(self):
        # n, A
        data = [(101, 5),
                (101, 3),
                (101, 1),
                (101, 0.05),
                (101, -0.05),
                (101, -1),
                (101, -3),
                (101, -5)]

        plt.figure()
        for i in range(len(data)):
            r = single_exponential(data[i][0], data[i][1])
            plt.plot(np.linspace(0.0, 1.0, data[i][0]), r, label='A={}'.format(data[i][1]))
        plt.xlabel('u')
        plt.ylabel('f(u)')
        plt.title('Single Exponential')
        plt.legend()
        plt.show()
        self.assertTrue(True)

    def test_double_exponential(self):
        # n, A1, A2, A3, A2取负值使得中间较密, 取正值使得中间稀疏，两边密集
        data = [(401, 0.5, -1.5, 0.5),
                (401, 0.5, 1.5, 0.5)]

        plt.figure()
        for i in range(len(data)):
            r = double_exponential(data[i][0], data[i][1], data[i][2], data[i][3])
            plt.plot(np.linspace(0.0, 1.0, data[i][0]), r, label='A1={}, A2={}, A3={}'.format(data[i][1], data[i][2], data[i][3]))
        plt.xlabel('u')
        plt.ylabel('f(u)')
        plt.title('Double Exponential')
        plt.legend()
        plt.show()
        self.assertTrue(True)

    def test_hyperbolic_tangent(self):
        # n, B
        data = [(201, 3), (201, 2), (201, 1)]

        plt.figure()
        for k, dt in enumerate(data):
            r = hyperbolic_tangent(dt[0], dt[1])
            plt.plot(np.linspace(0.0, 1.0, dt[0]), r, label='B={}'.format(dt[1]))
        plt.xlabel('u')
        plt.ylabel('f(u)')
        plt.title('Hyperbolic Tangent')
        plt.legend()
        plt.show()
        self.assertTrue(True)

    def test_hyperbolic_sine(self):
        # n, C
        data = [(201, 3), (201, 2), (201, 1)]

        plt.figure()
        for k, dt in enumerate(data):
            r = hyperbolic_sine(dt[0], dt[1])
            plt.plot(np.linspace(0.0, 1.0, dt[0]), r, label='C={}'.format(dt[1]))
        plt.xlabel('u')
        plt.ylabel('f(u)')
        plt.title('Hyperbolic Sine')
        plt.legend()
        plt.show()
        self.assertTrue(True)

    def test_newton_raphson(self):
        # pa1, pa2, pa3
        data = [(0.4, -1.2, 0.5)]
        ans = [0]

        self.assertTrue(len(data) == len(ans))

        for k, dt in enumerate(data):
            pa1, pa2, pa3 = dt
            p = (1 - pa1) * pa3 * (scipy.exp(pa2) - 1) / ((1 - pa3) * pa1 * pa2 * scipy.exp(pa2))
            p1z = math.log(p)

            def f(x):
                return scipy.exp(x) - 1.0 - p * x

            def pf(x):
                return scipy.exp(x) - p

            def ppf(x):
                return scipy.exp(x)

            fp1z = f(p1z)
            pa4 = optimize.newton(f, 5 * p1z, fprime=pf, maxiter=20, fprime2=ppf)
            print(p, p1z, fp1z, pa4, f(pa4))


if __name__ == '__main__':
    unittest.main()
