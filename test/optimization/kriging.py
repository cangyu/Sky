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

    def show(self, x1_rg=(-5, 10), x2_rg=(0, 15), x1_num=100, x2_num=100):
        """
        Show the function within specified range.
        :param x1_rg: Range of 'x1' coordinate.
        :param x2_rg: Range of 'x2' coordinate.
        :param x1_num: Number of points in 'x1' direction.
        :param x2_num: Number of points in 'x2' direction.
        :return: None.
        """

        '''Make data'''
        x1 = linear_expand(uniform(x1_num), x1_rg[0], x1_rg[1])
        x2 = linear_expand(uniform(x2_num), x2_rg[0], x2_rg[1])
        X1, X2 = np.meshgrid(x1, x2, indexing='ij')
        Y = np.empty((x1_num, x2_num), float)
        for i in range(x1_num):
            for j in range(x2_num):
                Y[i][j] = self.__call__(x1[i], x2[j])

        '''Build surf'''
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(X1, X2, Y, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        ax.zaxis.set_major_locator(LinearLocator(10))

        '''Add a color bar'''
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()


branin_func = Branin()
branin_func.show()

x = np.array([[0, 1, 2, 3, 4],
              [0, 1, 2, 3, 4]])

y = branin_func(x[0], x[1])

print(y)
