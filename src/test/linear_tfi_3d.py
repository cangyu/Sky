import numpy as np
import unittest
import math
from src.msh.linear_tfi import Linear_TFI_3D


def show(msh, pu, pv, pw):
    grid = msh.calc_msh(pu, pv, pw)
    # Linear_TFI_3D.show_msh(grid)
    print(grid)


def cuboid(L, W, H, I, J, K):
    """
    长方体网格
    :param L: 长 
    :param W: 宽
    :param H: 高
    :param I: L方向分段数
    :param J: W方向分段数
    :param K: H方向分段数
    :return: None
    """

    u_list = np.linspace(0, 1, I + 1)
    v_list = np.linspace(0, 1, J + 1)
    w_list = np.linspace(0, 1, K + 1)

    s1 = lambda v, w: np.array([0, v * W, w * H])
    s2 = lambda v, w: np.array([L, v * W, w * H])
    s3 = lambda w, u: np.array([u * L, 0, w * H])
    s4 = lambda w, u: np.array([u * L, W, w * H])
    s5 = lambda u, v: np.array([u * L, v * W, 0])
    s6 = lambda u, v: np.array([u * L, v * W, H])
    tmsh = Linear_TFI_3D(s1, s2, s3, s4, s5, s6)

    show(tmsh, u_list, v_list, w_list)


class T3L_Test(unittest.TestCase):
    """
    3D Linear Transfinite interpolation Test
    """

    def test_cuboid(self):
        cuboid(5, 4, 3, 10, 10, 10)
