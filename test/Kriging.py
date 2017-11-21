import unittest
import math
import numpy as np
from src.kriging import LHS, Kriging


class Branin(object):
    def __init__(self, a=1.0, b=5.1 / (4 * np.pi ** 2), c=5 / math.pi, r=6.0, s=10.0, t=1 / (8 * math.pi)):
        """
        Branin-Hoo function which has 3 global minima.
        """
        self.a = a
        self.b = b
        self.c = c
        self.r = r
        self.s = s
        self.t = t

    def __call__(self, x):
        """
        Calculate the value at (x1, x2).
        :param x: Parameter vector.
        :return: Value at (x1, x2).
        :rtype: float
        """

        x1 = x[0]
        x2 = x[1]
        ans = self.s
        ans += self.s * (1.0 - self.t) * math.cos(x1)
        ans += self.a * (x2 - self.b * x1 ** 2 + self.c * x1 - self.r) ** 2
        return ans


class KrigingTestCase(unittest.TestCase):
    def test_interp(self):
        branin_func = Branin()
        lhc = LHS(np.array([np.linspace(-5, 10, 20), np.linspace(0, 15, 20)]))

        x = lhc.sample()
        y = np.empty(len(x), float)
        for k, vx in enumerate(x):
            y[k] = branin_func(vx)

        kg = Kriging(x, y)
        x0 = (-math.pi, 12.275)
        x1 = (math.pi, 2.275)
        x2 = (9.42478, 2.475)
        print("\nGlobal minimum: {} Kriging: {} Actual: {}".format(x0, kg.interp(x0), branin_func(x0)))
        print("\nGlobal minimum: {} Kriging: {} Actual: {}".format(x1, kg.interp(x1), branin_func(x1)))
        print("\nGlobal minimum: {} Kriging: {} Actual: {}".format(x2, kg.interp(x2), branin_func(x2)))
        self.assertTrue(np.allclose(list(map(lambda u: kg.interp(u), x)), y))


if __name__ == '__main__':
    unittest.main()
