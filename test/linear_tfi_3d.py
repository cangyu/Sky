import numpy as np
import unittest
import math
from src.msh.linear_tfi import Linear_TFI_3D


def write_uniform_p3d(msh: Linear_TFI_3D, U: int, V: int, W: int, fn="msh_p3d.xyz"):
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    w_list = np.linspace(0, 1.0, W + 1)
    msh.calc_msh(u_list, v_list, w_list)
    msh.write_plot3d(fn)


def cuboid(L, W, H):
    """
    长方体网格
    :param L: 长 
    :param W: 宽
    :param H: 高
    :return: TFI_3D Mesh.
    """

    s1 = lambda v, w: np.array([0, v * W, w * H])
    s2 = lambda v, w: np.array([L, v * W, w * H])
    s3 = lambda w, u: np.array([u * L, 0, w * H])
    s4 = lambda w, u: np.array([u * L, W, w * H])
    s5 = lambda u, v: np.array([u * L, v * W, 0])
    s6 = lambda u, v: np.array([u * L, v * W, H])

    return Linear_TFI_3D(s1, s2, s3, s4, s5, s6)


class T3L_WritePlot3D_Test(unittest.TestCase):
    """
    3D Linear Transfinite Interpolation Test
    """

    def test_cuboid(self):
        U, V, W = 10, 8, 4
        msh = cuboid(5, 4, 3)
        write_uniform_p3d(msh, U, V, W, "cuboid.xyz")
