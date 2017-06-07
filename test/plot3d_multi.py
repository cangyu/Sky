import numpy as np
import math
from src.aircraft.wing import Airfoil
from src.nurbs.spline import Spline
from src.msh.tfi import Linear_TFI_2D
from src.msh.elliptic import Laplace2D
from src.msh.plot3d import PLOT3D, PLOT3D_Block


def airfoil(foil, L, R):
    """
    计算翼型流场边界
    :param foil: 翼型名称 
    :param L: 拉伸长度
    :param R: 远场半径
    :return: Parametric descriptions of mesh boundary.
    """

    '''Read airfoil data'''
    af = Airfoil(foil)
    for i in range(len(af.pts)):
        af.pts[i] *= L

    '''Parametric representation of airfoil'''
    crv = Spline()
    crv.interpolate(af.pts)
    cx = crv.x
    cy = crv.y
    c1 = lambda u: np.array([cx(u), cy(u), 0])

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

    '''Farfiled boundary'''
    sp = c2(1.0)
    ep = c4(1.0)
    sa = math.atan2(sp[1], sp[0])
    ea = math.atan2(ep[1], ep[0])
    if ea < 0:
        ea += math.pi * 2

    c3 = lambda u: np.array([r * math.cos((1 - u) * sa + u * ea), r * math.sin((1 - u) * sa + u * ea), 0])

    '''Tail'''
    c5 = lambda u: np.array([(1 - u) * af.pts[-1][0] + u * af.pts[0][0], (1 - u) * af.pts[-1][1] + u * af.pts[0][1], 0])
    ea2 = ea - math.pi * 2 if ea > 0 else ea
    c6 = lambda u: np.array([r * math.cos((1 - u) * ea2 + u * sa), r * math.sin((1 - u) * ea2 + u * sa), 0])

    return c1, c2, c3, c4, c5, c6


def write_uniform_airfoil(foil, L, R, U1, U2, V, fn="", delta_zeta=1.0, delta_eta=1.0):
    u1_list = np.linspace(0, 1.0, U1 + 1)
    u2_list = np.linspace(0, 1.0, U2 + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    c1, c2, c3, c4, c5, c6 = airfoil(foil, L, R)

    grid1 = Linear_TFI_2D(c1, c2, c3, c4)
    ppu, ppv = np.meshgrid(u1_list, v_list, indexing='ij')
    grid1.calc_grid(ppu, ppv)
    grid1 = Laplace2D(grid1.get_grid())
    grid1.calc_grid()

    grid2 = Linear_TFI_2D(c2, c5, c4, c6)
    ppu, ppv = np.meshgrid(u2_list, v_list, indexing='ij')
    grid2.calc_grid(ppu, ppv)

    if fn == "":
        fn += foil
        fn += "_{}_{}_{}_{}_{}_{}_{}_Multi.xyz".format(L, R, U1, U2, V, delta_zeta, delta_eta)

    msh = PLOT3D()
    msh.add_block(PLOT3D_Block.build_from_2d(grid1.get_grid()))
    msh.add_block(PLOT3D_Block.build_from_2d(grid2.get_grid()))
    msh.write(fn)


"""
Plot3D 多块网格测试
"""
if __name__ == '__main__':
    write_uniform_airfoil("NACA0012", 15, 50, 60, 10, 30)
