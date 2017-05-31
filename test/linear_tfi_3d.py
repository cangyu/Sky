import numpy as np
import unittest
import math
from src.msh.tfi import Linear_TFI_3D


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
    :return: TFI_3D Mesh
    """

    s1 = lambda v, w: np.array([0, v * W, w * H])
    s2 = lambda v, w: np.array([L, v * W, w * H])
    s3 = lambda w, u: np.array([u * L, 0, w * H])
    s4 = lambda w, u: np.array([u * L, W, w * H])
    s5 = lambda u, v: np.array([u * L, v * W, 0])
    s6 = lambda u, v: np.array([u * L, v * W, H])

    return Linear_TFI_3D(s1, s2, s3, s4, s5, s6)


def sect(R_MIN, R_MAX, THETA_MIN, THETA_MAX, H_MIN, H_MAX):
    """
    扇柱
    :param R_MIN: 内径
    :param R_MAX: 外径
    :param THETA_MIN: 起始角 
    :param THETA_MAX: 终止角
    :param H_MIN: 底部高度
    :param H_MAX: 顶部高度
    :return: TFI_3D Mesh
    """

    share = lambda p, a, b: (1 - p) * a + p * b

    s1 = lambda v, w: np.array([R_MIN * math.cos(math.radians(share(v, THETA_MIN, THETA_MAX))),
                                R_MIN * math.sin(math.radians(share(v, THETA_MIN, THETA_MAX))),
                                share(w, H_MIN, H_MAX)])

    s2 = lambda v, w: np.array([R_MAX * math.cos(math.radians(share(v, THETA_MIN, THETA_MAX))),
                                R_MAX * math.sin(math.radians(share(v, THETA_MIN, THETA_MAX))),
                                (1 - w) * H_MIN + w * H_MAX])

    s3 = lambda w, u: np.array([share(u, R_MIN, R_MAX) * math.cos(math.radians(THETA_MIN)),
                                share(u, R_MIN, R_MAX) * math.sin(math.radians(THETA_MIN)),
                                share(w, H_MIN, H_MAX)])

    s4 = lambda w, u: np.array([share(u, R_MIN, R_MAX) * math.cos(math.radians(THETA_MAX)),
                                share(u, R_MIN, R_MAX) * math.sin(math.radians(THETA_MAX)),
                                share(w, H_MIN, H_MAX)])

    s5 = lambda u, v: np.array([share(u, R_MIN, R_MAX) * math.cos(math.radians(share(v, THETA_MIN, THETA_MAX))),
                                share(u, R_MIN, R_MAX) * math.sin(math.radians(share(v, THETA_MIN, THETA_MAX))),
                                H_MIN])

    s6 = lambda u, v: np.array([share(u, R_MIN, R_MAX) * math.cos(math.radians(share(v, THETA_MIN, THETA_MAX))),
                                share(u, R_MIN, R_MAX) * math.sin(math.radians(share(v, THETA_MIN, THETA_MAX))),
                                H_MAX])

    return Linear_TFI_3D(s1, s2, s3, s4, s5, s6)


class T3L_WritePlot3D_Test(unittest.TestCase):
    """
    3D Linear Transfinite Interpolation Test
    """

    def test_cuboid(self):
        U, V, W = 10, 8, 4
        msh = cuboid(5, 4, 3)
        write_uniform_p3d(msh, U, V, W, "cuboid.xyz")

    def test_sect(self):
        U, V, W = 30, 15, 60
        msh = sect(5, 20, 60, 120, -50, 50)
        write_uniform_p3d(msh, U, V, W, "sect.xyz")
