import unittest
import numpy as np
import matplotlib.pyplot as plt
from src.msh.spacing import single_exponential, double_exponential


def plot_single(begin, end, N, A):
    rho = np.linspace(begin, end, N)
    r = single_exponential(begin, end, N, A)
    plt.figure()
    plt.plot(rho, r)
    plt.show()


def plot_double(begin, end, N, A1, A2, A3):
    rho = np.linspace(begin, end, N)
    r = double_exponential(begin, end, N, A1, A2, A3)
    plt.figure()
    plt.plot(rho, r)
    plt.show()


class SpacingTest(unittest.TestCase):
    @classmethod
    def test_single(cls):
        plot_single(0.0, 1.0, 101, 2)

    @classmethod
    def test_double(cls):
        plot_double(0.0, 1.0, 401, 0.4, -1.8, 0.46)
