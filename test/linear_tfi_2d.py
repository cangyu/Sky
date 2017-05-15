import unittest
import numpy as np
import math
import pylab
from src.msh.linear_tfi import Linear_TFI_2D
from src.aircraft.wing import Airfoil
from src.nurbs.curve import Spline
import matplotlib.pyplot as plt


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


def airfoil(foil, L, R):
    """
    翼型单块网格
    :param foil: 翼型名称 
    :param L: 拉伸长度
    :param R: 远场半径
    :return: TFI_2D mesh
    """

    '''Read airfoil data'''
    af = Airfoil(foil)
    pts = af.get_pts()

    '''Extend'''
    for p in pts:
        p *= L

    '''Parametric representation of airfoil'''
    crv = Spline(pts)
    cx = crv.x_rep()
    cy = crv.y_rep()
    c1 = lambda u: np.array([cx(u), cy(u)])

    '''Vertical boundary'''
    calc_dir = lambda start, end: (end[0] - start[0], end[1] - start[1])
    calc_len = lambda dir: math.pow(math.pow(dir[0], 2) + math.pow(dir[1], 2), 0.5)
    inner_product = lambda d1, d2: d1[0] * d2[0] + d1[1] * d2[1]

    dir1 = calc_dir(pts[1], pts[0])
    len1 = calc_len(dir1)
    dir2 = calc_dir(pts[-2], pts[-1])
    len2 = calc_len(dir2)

    rotate = lambda dir, theta: (dir[0] * math.cos(math.radians(theta)) - dir[1] * math.sin(math.radians(theta)),
                                 dir[0] * math.sin(math.radians(theta)) + dir[1] * math.cos(math.radians(theta)))
    dir1 = rotate(dir1, 90)
    dir2 = rotate(dir2, 270)
    dir3 = (dir1[0] / len1, dir1[1] / len1)
    dir4 = (dir2[0] / len2, dir2[1] / len2)

    c2 = lambda v: np.array([pts[0][0] + v * R * dir3[0], pts[0][1] + v * R * dir3[1]])
    r = calc_len(c2(1.0))

    dir5 = (pts[-1][0], pts[-1][1])
    l4 = calc_len(dir5)
    theta = math.pi - math.acos(inner_product(dir4, dir5) / l4)
    alpha = math.asin(l4 / r * math.sin(theta))
    beta = math.pi - theta - alpha
    b = r / math.sin(theta) * math.sin(beta)

    c4 = lambda v: np.array([pts[-1][0] + v * b * dir4[0], pts[-1][1] + v * b * dir4[1]])

    '''Farfiled boundary'''
    sp = c2(1.0)
    ep = c4(1.0)
    sa = math.atan2(sp[1], sp[0])
    ea = math.atan2(ep[1], ep[0])
    if ea < 0:
        ea += math.pi * 2

    c3 = lambda u: np.array([r * math.cos((1 - u) * sa + u * ea), r * math.sin((1 - u) * sa + u * ea)])

    return Linear_TFI_2D(c1, c2, c3, c4)


class T2L_Show_Test(unittest.TestCase):
    """
    Test Python representation.
    """

    def test_rectangular(self):
        msh = rectangular(5, 4)
        U, V = 10, 8
        show_uniform(msh, U, V)

    def test_circle(self):
        msh = circle(1, 2)
        U, V = 5, 10
        show_uniform(msh, U, V)

    def test_eccentic(self):
        msh = eccentric_circle(-10, 4, 25)
        U, V = 15, 40
        show_uniform(msh, U, V)

    def test_crv_rect(self):
        msh = curve_rect(100, 40, 60, 10)
        U, V = 50, 25
        show_uniform(msh, U, V)

    def test_airfoil(self):
        msh = airfoil("NACA0012", 10, 50)
        U, V = 60, 15
        show_uniform(msh, U, V)


class T2L_WritePlot3D_Test(unittest.TestCase):
    """
    Test Plot3D file generation.
    """

    def test_rectangular(self):
        msh = rectangular(5, 4)
        U, V = 10, 8
        write_uniform_p3d(msh, U, V, "rect.xyz")

    def test_circle(self):
        msh = circle(1, 2)
        U, V = 5, 10
        write_uniform_p3d(msh, U, V, "circle.xyz")

    def test_eccentic(self):
        msh = eccentric_circle(-10, 4, 25)
        U, V = 15, 40
        write_uniform_p3d(msh, U, V, "eccentric.xyz")

    def test_crv_rect(self):
        msh = curve_rect(100, 40, 60, 10)
        U, V = 50, 25
        write_uniform_p3d(msh, U, V, "crv_rect.xyz")

    def test_airfoil(self):
        msh = airfoil("NACA0012", 10, 50)
        U, V = 60, 15
        write_uniform_p3d(msh, U, V, "NACA0012.xyz")
