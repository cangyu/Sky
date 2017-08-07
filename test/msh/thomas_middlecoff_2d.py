import numpy as np
import math
from src.aircraft.wing import Airfoil
from src.nurbs.curve import GlobalInterpolatedCrv
from src.msh.elliptic import ThomasMiddlecoff2D
from src.msh.plot3d import PLOT3D_Block, PLOT3D
from src.msh.tfi import LinearTFI2D
from src.msh.spacing import single_exponential, double_exponential


def write_airfoil_o_msh(foil, L, R, U, V, fn=''):
    """
    计算翼型O型流场边界
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

    '''farfield boundary'''
    sp = c2(1.0)
    ep = c4(1.0)
    sa = math.atan2(sp[1], sp[0])
    ea = math.atan2(ep[1], ep[0])
    if ea < 0:
        ea += math.pi * 2

    c3 = lambda u: np.array([r * math.cos((1 - u) * sa + u * ea), r * math.sin((1 - u) * sa + u * ea), 0])

    if fn == '':
        fn = foil + "_{}_{}_{}_{}_{}_{}_TM.xyz".format(L, R, U, V, 1.0, 1.0)

    tfi_grid = LinearTFI2D(c1, c2, c3, c4)

    u_list = double_exponential(U + 1, 0.5, -1.5, 0.5)
    v_list = single_exponential(V + 1, 1.5)
    tfi_grid.calc_grid(u_list, v_list)

    tm_grid = ThomasMiddlecoff2D(tfi_grid.get_grid())
    tm_grid.calc_grid()
    msh = PLOT3D()
    msh.add_block(PLOT3D_Block.build_from_2d(tm_grid.get_grid()))
    msh.write(fn)


if __name__ == '__main__':
    U, V = 400, 125
    write_airfoil_o_msh("NACA0012", 10, 50, U, V)
