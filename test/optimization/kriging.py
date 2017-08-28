import numpy as np
from src.msh.spacing import uniform, linear_expand
from src.opt.surrogate import Kriging
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator


class Branin(object):
    def __init__(self, a=1.0, b=5.1 / (4 * np.pi ** 2), c=5 / np.pi, r=6.0, s=10.0, t=1 / (8 * np.pi)):
        """
        Branin-Hoo function which has 3 global minima.
        """
        self.a = a
        self.b = b
        self.c = c
        self.r = r
        self.s = s
        self.t = t

    def __call__(self, x1, x2):
        """
        Calculate the value at (x1, x2).
        :param x1: First coordinate component.
        :type x1: float
        :param x2: Second coordinate component.
        :type x2: float
        :return: Value at (x1, x2).
        :rtype: float
        """

        ans = self.s
        ans += self.s * (1.0 - self.t) * np.cos(x1)
        ans += self.a * (x2 - self.b * x1 ** 2 + self.c * x1 - self.r) ** 2
        return ans


branin_func = Branin()

x = np.array([[0, 1, 2, 3, 4],
              [0, 1, 2, 3, 4]])

y = branin_func(x[0], x[1])

print(x.transpose())
print(y)

kg = Kriging(x.transpose(), y)
