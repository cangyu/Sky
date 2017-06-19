import unittest
import numpy as np
import matplotlib.pyplot as plt
from src.msh.spacing import single_exponential, double_exponential, hyperbolic_tangent, hyperbolic_sine


def plot_single(N, A):
    r = single_exponential(N, A)
    plt.figure()
    plt.plot(np.linspace(0.0, 1.0, N), r)
    plt.xlabel('u')
    plt.ylabel('f(u)')
    plt.title('Single Exponential\nA={}'.format(A))
    plt.show()


def plot_double(N, A1, A2, A3):
    r = double_exponential(N, A1, A2, A3)
    plt.figure()
    plt.plot(np.linspace(0.0, 1.0, N), r)
    plt.xlabel('u')
    plt.ylabel('f(u)')
    plt.title('Double Exponential\nA1={}, A2={}, A3={}'.format(A1, A2, A3))
    plt.show()


def plot_hyperbolic_tangent(N, B):
    r = hyperbolic_tangent(N, B)
    plt.figure()
    plt.plot(np.linspace(0.0, 1.0, N), r)
    plt.xlabel('u')
    plt.ylabel('f(u)')
    plt.title('Hyperbolic Tangent\nB={}'.format(B))
    plt.show()


def plot_hyperbolic_sine(N, C):
    r = hyperbolic_sine(N, C)
    plt.figure()
    plt.plot(np.linspace(0.0, 1.0, N), r)
    plt.xlabel('u')
    plt.ylabel('f(u)')
    plt.title('Hyperbolic Tangent\nC={}'.format(C))
    plt.show()


class SpacingTest(unittest.TestCase):
    @classmethod
    def test(cls):
        plot_single(101, 2)
        plot_double(401, 0.5, -1.5, 0.5)  # A2取负值使得中间较密
        plot_double(401, 0.5, 1.5, 0.5)  # A2取正值使得中间稀疏，两边密集
        plot_hyperbolic_tangent(201, 2)
        plot_hyperbolic_sine(201, 2.5)
