import numpy as np
import unittest
import math
from src.msh.tfi import LinearTFI3D
from src.msh.fluent import XF_MSH


def sect(r_min, r_max, theta_min, theta_max, h_min, h_max, nu, nv, nw):
    """
    扇柱单块网格
    :param r_min: 内径
    :type r_min: float
    :param r_max: 外径
    :type r_max: float
    :param theta_min: 起始角
    :type theta_min: float
    :param theta_max: 终止角
    :type theta_max: float
    :param h_min: 底部高度
    :type h_min: float
    :param h_max: 顶部高度
    :type h_max: float
    :param nu: X方向网格数
    :type nu: int
    :param nv: Y方向网格数
    :type nv: int
    :param nw: Z方向网格数
    :type nw: int
    :return: TFI_3D Mesh
    """

    def share(p, a, b):
        return (1 - p) * a + p * b

    def s1(v, w):
        return np.array([r_min * math.cos(math.radians(share(v, theta_min, theta_max))),
                         r_min * math.sin(math.radians(share(v, theta_min, theta_max))),
                         share(w, h_min, h_max)])

    def s2(v, w):
        return np.array([r_max * math.cos(math.radians(share(v, theta_min, theta_max))),
                         r_max * math.sin(math.radians(share(v, theta_min, theta_max))),
                         (1 - w) * h_min + w * h_max])

    def s3(w, u):
        return np.array([share(u, r_min, r_max) * math.cos(math.radians(theta_min)),
                         share(u, r_min, r_max) * math.sin(math.radians(theta_min)),
                         share(w, h_min, h_max)])

    def s4(w, u):
        return np.array([share(u, r_min, r_max) * math.cos(math.radians(theta_max)),
                         share(u, r_min, r_max) * math.sin(math.radians(theta_max)),
                         share(w, h_min, h_max)])

    def s5(u, v):
        return np.array([share(u, r_min, r_max) * math.cos(math.radians(share(v, theta_min, theta_max))),
                         share(u, r_min, r_max) * math.sin(math.radians(share(v, theta_min, theta_max))),
                         h_min])

    def s6(u, v):
        return np.array([share(u, r_min, r_max) * math.cos(math.radians(share(v, theta_min, theta_max))),
                         share(u, r_min, r_max) * math.sin(math.radians(share(v, theta_min, theta_max))),
                         h_max])

    msh = LinearTFI3D(s1, s2, s3, s4, s5, s6)

    u_list = np.linspace(0, 1.0, nu + 1)
    v_list = np.linspace(0, 1.0, nv + 1)
    w_list = np.linspace(0, 1.0, nw + 1)
    msh.calc_grid(u_list, v_list, w_list)

    return msh.get_grid()


class TestFromStr3D(unittest.TestCase):
    @staticmethod
    def test():
        grid = sect(10, 32, 0, 160, -15, 15, 25, 30, 20)
        fluent_msh = XF_MSH.from_str3d(grid)
        fluent_msh.save("sect.msh")
