import unittest
import numpy as np
import math
from src.msh.tfi import LinearTFI2D
from src.msh.plot3d import PLOT3D_Block, PLOT3D
from src.aircraft.wing import Airfoil
from src.geom.curve import GlobalInterpolatedCrv


def write_uniform_p3d(msh: LinearTFI2D, U: int, V: int, fn="msh_p3d.xyz"):
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    msh.calc_grid(u_list, v_list)
    grid = PLOT3D()
    grid.add_block(PLOT3D_Block.build_from_2d(msh.get_grid()))
    grid.write(fn)


def rectangular(L: float, W: float):
    c1 = lambda u: np.array([L * u, 0, 0])
    c3 = lambda u: np.array([L * u, W, 0])
    c2 = lambda v: np.array([0, W * v, 0])
    c4 = lambda v: np.array([L, W * v, 0])
    return LinearTFI2D(c1, c2, c3, c4)


def circle(R1: float, R2: float):
    c1 = lambda u: np.array([(1 - u) * R1 + u * R2, 0, 0])
    c3 = lambda u: np.array([0, (1 - u) * R1 + u * R2, 0])
    c2 = lambda v: np.array([R1 * math.cos(0.5 * math.pi * v), R1 * math.sin(0.5 * math.pi * v), 0])
    c4 = lambda v: np.array([R2 * math.cos(0.5 * math.pi * v), R2 * math.sin(0.5 * math.pi * v), 0])
    return LinearTFI2D(c1, c2, c3, c4)


def eccentric_circle(delta: float, R1: float, R2: float):
    c1 = lambda u: np.array([(delta + R1) * (1 - u) + R2 * u, 0, 0])
    c3 = lambda u: np.array([(delta - R1) * (1 - u) - R2 * u, 0, 0])
    c2 = lambda v: np.array([R1 * math.cos(math.pi * v) + delta, R1 * math.sin(math.pi * v), 0])
    c4 = lambda v: np.array([R2 * math.cos(math.pi * v), R2 * math.sin(math.pi * v), 0])
    return LinearTFI2D(c1, c2, c3, c4)


def curve_rect(L: float, H1: float, H2: float, H3: float):
    c1 = lambda u: np.array([u * L, 4 * H3 * u * (1 - u), 0])
    c3 = lambda u: np.array([u * L, (H1 * (1 - u * u) + H2 * u * u), 0])
    c2 = lambda v: np.array([0, v * H1, 0])
    c4 = lambda v: np.array([L, v * H2, 0])
    return LinearTFI2D(c1, c2, c3, c4)


def airfoil(foil, L, R):
    """
    翼型单块网格
    :param foil: 翼型名称 
    :param L: 拉伸长度
    :param R: 远场半径
    :return: TFI_2D mesh
    """

    '''Read airfoil data'''
    af = Airfoil.from_database(foil)

    '''Extend'''
    for i in range(len(af.pts)):
        af.pts[i] *= L

    '''Parametric representation of airfoil'''
    crv = GlobalInterpolatedCrv(af.pts, 5)
    c1 = lambda u: crv(u)

    '''Vertical boundary'''
    calc_dir = lambda start, end: (end[0] - start[0], end[1] - start[1])
    calc_len = lambda dir: math.pow(math.pow(dir[0], 2) + math.pow(dir[1], 2), 0.5)
    inner_product = lambda d1, d2: d1[0] * d2[0] + d1[1] * d2[1]

    dir1 = calc_dir(af.pts[1], af.pts[0])
    len1 = calc_len(dir1)
    dir2 = calc_dir(af.pts[-2], af.pts[-1])
    len2 = calc_len(dir2)

    rotate = lambda dir, theta: (dir[0] * math.cos(math.radians(theta)) - dir[1] * math.sin(math.radians(theta)),
                                 dir[0] * math.sin(math.radians(theta)) + dir[1] * math.cos(math.radians(theta)))
    dir1 = rotate(dir1, 90)
    dir2 = rotate(dir2, 270)
    dir3 = (dir1[0] / len1, dir1[1] / len1)
    dir4 = (dir2[0] / len2, dir2[1] / len2)

    c2 = lambda v: np.array([af.pts[0][0] + v * R * dir3[0], af.pts[0][1] + v * R * dir3[1], 0])
    r = calc_len(c2(1.0))

    dir5 = (af.pts[-1][0], af.pts[-1][1])
    l4 = calc_len(dir5)
    theta = math.pi - math.acos(inner_product(dir4, dir5) / l4)
    alpha = math.asin(l4 / r * math.sin(theta))
    beta = math.pi - theta - alpha
    b = r / math.sin(theta) * math.sin(beta)

    c4 = lambda v: np.array([af.pts[-1][0] + v * b * dir4[0], af.pts[-1][1] + v * b * dir4[1], 0])

    '''Farfield boundary'''
    sp = c2(1.0)
    ep = c4(1.0)
    sa = math.atan2(sp[1], sp[0])
    ea = math.atan2(ep[1], ep[0])
    if ea < 0:
        ea += math.pi * 2

    c3 = lambda u: np.array([r * math.cos((1 - u) * sa + u * ea), r * math.sin((1 - u) * sa + u * ea), 0])

    return LinearTFI2D(c1, c2, c3, c4)


class T2L_WritePlot3D_Test(unittest.TestCase):
    """
    Test Plot3D file generation.
    """

    @staticmethod
    def test_rectangular():
        msh = rectangular(5, 4)
        U, V = 10, 8
        write_uniform_p3d(msh, U, V, "rect.xyz")

    @staticmethod
    def test_circle():
        msh = circle(1, 2)
        U, V = 5, 10
        write_uniform_p3d(msh, U, V, "circle.xyz")

    @staticmethod
    def test_eccentric():
        msh = eccentric_circle(-10, 4, 25)
        U, V = 15, 40
        write_uniform_p3d(msh, U, V, "eccentric.xyz")

    @staticmethod
    def test_crv_rect():
        msh = curve_rect(100, 40, 60, 10)
        U, V = 50, 25
        write_uniform_p3d(msh, U, V, "crv_rect.xyz")

    @staticmethod
    def test_airfoil():
        msh = airfoil("NACA0012", 10, 50)
        U, V = 60, 15
        write_uniform_p3d(msh, U, V, "NACA0012.xyz")
