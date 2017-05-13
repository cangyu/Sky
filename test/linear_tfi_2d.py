import unittest
import numpy as np
import math
import pylab
from src.msh.linear_tfi import Linear_TFI_2D


def show(msh: Linear_TFI_2D, pu, pv):
    grid = msh.calc_msh(pu, pv)
    U, V, D = grid.shape
    x = np.zeros((V, U))
    y = np.zeros((V, U))

    for i in range(0, V):
        for j in range(0, U):
            x[i][j] = grid[j][i][0]
            y[i][j] = grid[j][i][1]

    pylab.plot(x, y)
    pylab.plot(np.vstack((x[:, 0], x[:, -1])), np.vstack((y[:, 0], y[:, -1])))
    pylab.axis('scaled')
    pylab.show()


def show_uniform(msh: Linear_TFI_2D, U: int, V: int):
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    show(msh, u_list, v_list)


def write_uniform_p3d(msh: Linear_TFI_2D, U: int, V: int, fn="msh_p3d.xyz"):
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    msh.write_plot3d(u_list, v_list, fn)


def rectangular(L: float, W: float):
    c1 = lambda u: np.array([L * u, 0])
    c3 = lambda u: np.array([L * u, W])
    c2 = lambda v: np.array([0, W * v])
    c4 = lambda v: np.array([L, W * v])
    return Linear_TFI_2D(c1, c2, c3, c4)


def circle(R1: float, R2: float):
    c1 = lambda u: np.array([(1 - u) * R1 + u * R2, 0])
    c3 = lambda u: np.array([0, (1 - u) * R1 + u * R2])
    c2 = lambda v: np.array([R1 * math.cos(0.5 * math.pi * v), R1 * math.sin(0.5 * math.pi * v)])
    c4 = lambda v: np.array([R2 * math.cos(0.5 * math.pi * v), R2 * math.sin(0.5 * math.pi * v)])
    return Linear_TFI_2D(c1, c2, c3, c4)


def eccentric_circle(delta: float, R1: float, R2: float):
    c1 = lambda u: np.array([(delta + R1) * (1 - u) + R2 * u, 0])
    c3 = lambda u: np.array([(delta - R1) * (1 - u) - R2 * u, 0])
    c2 = lambda v: np.array([R1 * math.cos(math.pi * v) + delta, R1 * math.sin(math.pi * v)])
    c4 = lambda v: np.array([R2 * math.cos(math.pi * v), R2 * math.sin(math.pi * v)])
    return Linear_TFI_2D(c1, c2, c3, c4)


def curve_rect(L: float, H1: float, H2: float, H3: float):
    c1 = lambda u: np.array([u * L, 4 * H3 * u * (1 - u)])
    c3 = lambda u: np.array([u * L, (H1 * (1 - u * u) + H2 * u * u)])
    c2 = lambda v: np.array([0, v * H1])
    c4 = lambda v: np.array([L, v * H2])
    return Linear_TFI_2D(c1, c2, c3, c4)


class T2L_Test(unittest.TestCase):
    """
    2D Linear Transfinite Interpolation Test
    """

    def test_rectangular(self):
        msh = rectangular(5, 4)
        U, V = 10, 8
        show_uniform(msh, U, V)
        write_uniform_p3d(msh, U, V, "rect.xyz")

    def test_circle(self):
        msh = circle(1, 2)
        U, V = 5, 10
        show_uniform(msh, U, V)
        write_uniform_p3d(msh, U, V, "circle.xyz")

    def test_eccentic(self):
        msh = eccentric_circle(-10, 4, 25)
        U, V = 15, 40
        show_uniform(msh, U, V)
        write_uniform_p3d(msh, U, V, "eccentric.xyz")

    def test_crv_rect(self):
        msh = curve_rect(100, 40, 60, 10)
        U, V = 50, 25
        show_uniform(msh, U, V)
        write_uniform_p3d(msh, U, V, "crv_rect.xyz")
