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
    def test_single(self):
        plot_single(0.0, 1.0, 101, 2)

    def test_double(self):
        plot_double(0.0, 1.0, 401, 0.45, -1.2, 0.52)


if __name__ == '__main__':
    unittest.main()
