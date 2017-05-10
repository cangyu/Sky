import numpy as np
import unittest
import math
from src.msh.linear_tfi import Linear_TFI_2D


def show(msh, pu, pv):
    Linear_TFI_2D.show_msh(msh.calc_msh(pu, pv))


def rectangular(L, W, U, V):
    c1 = lambda u: np.array([L * u, 0])
    c3 = lambda u: np.array([L * u, W])
    c2 = lambda v: np.array([0, W * v])
    c4 = lambda v: np.array([L, W * v])

    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    tmsh = Linear_TFI_2D(c1, c2, c3, c4)

    show(tmsh, u_list, v_list)


def circle(R1, R2, U, V):
    c1 = lambda u: np.array([(1 - u) * R1 + u * R2, 0])
    c3 = lambda u: np.array([0, (1 - u) * R1 + u * R2])
    c2 = lambda v: np.array([R1 * math.cos(0.5 * math.pi * v), R1 * math.sin(0.5 * math.pi * v)])
    c4 = lambda v: np.array([R2 * math.cos(0.5 * math.pi * v), R2 * math.sin(0.5 * math.pi * v)])

    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    tmsh = Linear_TFI_2D(c1, c2, c3, c4)

    show(tmsh, u_list, v_list)


class T2L_Test(unittest.TestCase):
    """
    2D Linear Transfinite interpolation Test
    """

    def test_rectangular(self):
        rectangular(5, 4, 10, 8)

    def test_circle(self):
        circle(1, 2, 5, 10)
