import numpy as np
import unittest
import math
from src.msh.linear_tfi import Linear_TFI_3D


def show(msh, pu, pv, pw):
    grid = msh.calc_msh(pu, pv, pw)
    # Linear_TFI_3D.show_msh(grid)
    print(grid)

def show_msh(grid):
    U, V, W, D = grid.shape

    x = np.zeros((W, V, U))
    y = np.zeros((W, V, U))
    z = np.zeros((W, V, U))

    for k in range(0, W):
        for j in range(0, V):
            for i in range(0, U):
                x[k][j][i] = grid[i][j][k][0]
                y[k][j][i] = grid[i][j][k][1]
                z[k][j][i] = grid[i][j][k][2]

    pylab.plot(x, y, z)
    # pylab.plot(np.vstack((x[:, 0], x[:, -1])), np.vstack((y[:, 0], y[:, -1])))
    pylab.axis('scaled')
    pylab.show()



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
