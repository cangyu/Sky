import unittest
import numpy as np
from src.msh.tfi import LinearTFI2D
from src.msh.fluent import XF_MSH, BCType
from src.geom.curve import Arc, Line
from src.aircraft.wing import WingProfile


def rectangular(u, v, l, w):
    """
    矩形
    :param u: U方向网格数量
    :type u: int
    :param v: V方向网格数量
    :type v: int
    :param l: 矩形长度
    :type l: float
    :param w: 矩形宽度
    :type w: float
    :return: None
    """

    msh = LinearTFI2D(lambda u: np.array([l * u, 0, 0]),
                      lambda v: np.array([0, w * v, 0]),
                      lambda u: np.array([l * u, w, 0]),
                      lambda v: np.array([l, w * v, 0]))

    u_list = np.linspace(0, 1.0, u + 1)
    v_list = np.linspace(0, 1.0, v + 1)
    msh.calc_grid(u_list, v_list)

    fluent_grid = XF_MSH.from_str2d(msh.get_grid())
    fluent_grid.save('rect_{}_{}.msh'.format(u, v))


def curve_rect(U, V, L, H1, H2, H3):
    """
    曲边矩形
    :param U: U方向网格数量
    :type U: int
    :param V: V方向网格数量
    :type V: int
    :param L: 矩形长度
    :type L: float
    :param H1: 控制高度1
    :type H1: float
    :param H2: 控制高度2
    :type H2: float
    :param H3: 控制高度3
    :type H3: float
    :return: None
    """

    msh = LinearTFI2D(lambda u: np.array([u * L, 4 * H3 * u * (1 - u), 0]),
                      lambda v: np.array([0, v * H1, 0]),
                      lambda u: np.array([u * L, (H1 * (1 - u * u) + H2 * u * u), 0]),
                      lambda v: np.array([L, v * H2, 0]))

    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    msh.calc_grid(u_list, v_list)

    fluent_grid = XF_MSH.from_str2d(msh.get_grid())
    fluent_grid.save('crv_rect_{}_{}.msh'.format(U, V))


def airfoil(U, V, foil, ends, thk, order, arc_start, arc_end, theta, nv):
    """
    半圆翼型
    :param U: U方向网格数量
    :param V: V方向网格数量
    :param foil: 翼型名称
    :param ends: 翼型拉伸前后端点
    :param thk: 翼型厚度缩放因子
    :param order: 翼型曲线次数
    :param arc_start: 远场圆弧起始点
    :param arc_end: 远场圆弧终止点
    :param theta: 远场圆弧角度
    :param nv: 远场圆弧法向量
    :return: None
    """

    u0 = WingProfile(foil, ends, thk).nurbs_rep(order)
    u1 = Arc.from_2pnt(arc_start, arc_end, theta, nv)
    v0 = Line(u0.start, u1.start)
    v1 = Line(u0.end, u1.end)

    msh = LinearTFI2D(u0, v0, u1, v1)

    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    msh.calc_grid(u_list, v_list)

    fluent_grid = XF_MSH.from_str2d(msh.get_grid(), bc=(BCType.Wall, BCType.Outflow, BCType.VelocityInlet, BCType.Outflow))
    fn = '{}_{}_{}.msh'.format(foil, U, V)
    fluent_grid.save(fn)


class SingleFluentMeshTest(unittest.TestCase):
    @staticmethod
    def test():
        rectangular(30, 10, 3, 1)
        curve_rect(50, 25, 100, 40, 60, 10)
        airfoil(130, 45, 'M6', [[0, 0, 5], [20, 0, 5]], 1.6, 5, [25, 60, 5], [25, -60, 5], 180, [0, 0, 1])
