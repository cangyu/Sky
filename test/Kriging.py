import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.kriging import LHS, Kriging


class Branin(object):
    def __init__(self, a=1.0, b=5.1 / (4 * math.pi ** 2), c=5 / math.pi, r=6.0, s=10.0, t=1 / (8 * math.pi)):
        """
        Branin-Hoo function which has 3 global minimal.
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


def interp_compare():
    branin_func = Branin()
    sample_num = 20
    lhc = LHS(np.array([np.linspace(-5, 10, sample_num), np.linspace(0, 15, sample_num)]))

    x = lhc.sample()
    y = np.array([branin_func(vx) for vx in x])
    kg = Kriging(x, y)
    for _x in [[-math.pi, 12.275], [math.pi, 2.275], [9.42478, 2.475]]:
        print("\nGlobal minimum: {} Kriging: {} Actual: {}".format(_x, kg.interp(_x), branin_func(_x)))

    '''plot contours'''
    n = 200
    X1 = np.linspace(-5, 10, n)
    X2 = np.linspace(0, 15, n)
    X, Y = np.meshgrid(X1, X2)

    actual = np.empty((n, n), float)
    predict = np.empty((n, n), float)

    for i in range(n):
        for j in range(n):
            x = (X1[i], X2[j])
            actual[i][j] = branin_func(x)
            predict[i][j] = kg.interp(x)

    fig = plt.figure()
    actual_ax = fig.add_subplot(121)
    actual_ax.contourf(X, Y, actual, alpha=0.75)
    cs = actual_ax.contour(X, Y, actual, 12, colors='black', linewidth=0.2)
    actual_ax.clabel(cs, inline=True, fontsize=8)
    actual_ax.set_xlabel('x1')
    actual_ax.set_ylabel('x2')
    actual_ax.set_title('Branin-Hoo Function contour')

    predict_ax = fig.add_subplot(122)
    predict_ax.contourf(X, Y, predict, alpha=0.75)
    cs = predict_ax.contour(X, Y, predict, 12, colors='black', linewidth=0.2)
    predict_ax.clabel(cs, inline=True, fontsize=8)
    predict_ax.set_xlabel('x1')
    predict_ax.set_ylabel('x2')
    predict_ax.set_title('Kriging predictions contour')

    fig.tight_layout()
    fig.savefig('kriging_test.png', dpi=800)
    plt.show()


if __name__ == '__main__':
    interp_compare()
