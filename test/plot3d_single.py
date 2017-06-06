import unittest
import numpy as np
import math
from src.msh.plot3d import PLOT3D, PLOT3D_Block


def rect(X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX, U, V, W, fout=""):
    """
    生成一个长方体
    :param X_MIN: 长度方向最小值 
    :param X_MAX: 长度方向最大值
    :param Y_MIN: 宽度方向最小值
    :param Y_MAX: 宽度方向最大值
    :param Z_MIN: 高度方向最小值
    :param Z_MAX: 高度方向最大值
    :param U: 长度方向分段数
    :param V: 宽度方向分段数
    :param W: 高度方向分段数
    :param fout: 输出文件名
    :return: None
    """

    pts = np.zeros((U + 1, V + 1, W + 1, 3))

    u_list = np.linspace(X_MIN, X_MAX, U + 1)
    v_list = np.linspace(Y_MIN, Y_MAX, V + 1)
    w_list = np.linspace(Z_MIN, Z_MAX, W + 1)

    if fout == "":
        fout += "rect"
        fout += "_{}_{}".format(X_MIN, X_MAX)
        fout += "_{}_{}".format(Y_MIN, Y_MAX)
        fout += "_{}_{}".format(Z_MIN, Z_MAX)
        fout += "_{}_{}_{}.xyz".format(U, V, W)

    for i in range(0, U + 1):
        for j in range(0, V + 1):
            for k in range(0, W + 1):
                pts[i][j][k] = np.array([u_list[i], v_list[j], w_list[k]])

    msh = PLOT3D()
    msh.add_block(PLOT3D_Block.build_from_3d(pts))
    msh.write(fout)


def sect(R_MIN, R_MAX, THETA_MIN, THETA_MAX, H_MIN, H_MAX, U, V, W, fout=""):
    """
    生成一个扇形柱
    :param R_MIN: 扇形内径
    :param R_MAX: 扇形外径
    :param THETA_MIN: 起始角度(Degree)
    :param THETA_MAX: 终止角度(Degree)
    :param H_MIN: 最小高度
    :param H_MAX: 最大高度
    :param U: 径向分段数
    :param V: 周向分段数
    :param W: 高度方向分段数
    :param fout: 输出文件名
    :return: None
    """

    pts = np.zeros((U + 1, V + 1, W + 1, 3))

    u_list = np.linspace(R_MIN, R_MAX, U + 1)
    v_list = np.linspace(THETA_MIN, THETA_MAX, V + 1)
    w_list = np.linspace(H_MIN, H_MAX, W + 1)

    if fout == "":
        fout += "sect"
        fout += "_{}_{}".format(R_MIN, R_MAX)
        fout += "_{}_{}".format(THETA_MIN, THETA_MAX)
        fout += "_{}_{}".format(H_MIN, H_MAX)
        fout += "_{}_{}_{}.xyz".format(U, V, W)

    for i in range(0, U + 1):
        for j in range(0, V + 1):
            for k in range(0, W + 1):
                ct = math.radians(v_list[j])
                pts[i][j][k] = np.array([u_list[i] * math.cos(ct), u_list[i] * math.sin(ct), w_list[k]])

    msh = PLOT3D()
    msh.add_block(PLOT3D_Block.build_from_3d(pts))
    msh.write(fout)


class Plot3D_Test(unittest.TestCase):
    @staticmethod
    def test():
        rect(0, 100, 0, 60, 0, 40, 60, 15, 20)
        sect(50, 100, 60, 320, 0, 30, 60, 15, 20)
