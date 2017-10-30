import math

import numpy as np
from src.msh.plot3d import PLOT3D_Block, PLOT3D
from src.msh.tfi import LinearTFI2D

from nurbs import GlobalInterpolatedCrv
from spacing import single_exponential, double_exponential
from src.aircraft.wing import Airfoil
from src.msh.elliptic import Laplace2D


def write_airfoil_msh(foil, L, R, U, V, fn="", delta_zeta=1.0, delta_eta=1.0):
    """
    计算翼型流场边界
    :param foil: 翼型名称 
    :param L: 拉伸长度
    :param R: 远场半径
    :return: Parametric descriptions of mesh boundary.
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

    '''Farfiled boundary'''
    sp = c2(1.0)
    ep = c4(1.0)
    sa = math.atan2(sp[1], sp[0])
    ea = math.atan2(ep[1], ep[0])
    if ea < 0:
        ea += math.pi * 2

    c3 = lambda u: np.array([r * math.cos((1 - u) * sa + u * ea), r * math.sin((1 - u) * sa + u * ea), 0])

    u_list = double_exponential(U + 1, 0.5, -1.5, 0.5)
    v_list = single_exponential(V + 1, 3)
    grid = LinearTFI2D(c1, c2, c3, c4)
    grid.calc_grid(u_list, v_list)
    lgrid = Laplace2D(grid.get_grid())
    lgrid.calc_grid()

    if fn == "":
        fn += foil + "_{}_{}_{}_{}_{}_{}_Laplace.xyz".format(L, R, U, V, delta_zeta, delta_eta)

    msh = PLOT3D()
    msh.add_block(PLOT3D_Block.build_from_2d(lgrid.get_grid()))
    msh.write(fn)


"""
Test grid generation.
"""
if __name__ == '__main__':
    write_airfoil_msh('NACA0012', 10, 50, 60, 25)
    write_airfoil_msh('M6', 10, 50, 60, 25)
