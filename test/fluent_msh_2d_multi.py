import unittest
import numpy as np
import math
from src.nurbs.curve import GlobalInterpolatedCrv, Line, Arc
from src.aircraft.wing import WingProfile
from src.msh.spacing import single_exponential, double_exponential, hyperbolic_tangent
from src.msh.elliptic import Laplace2D, ThomasMiddlecoff2D
from src.msh.tfi import Linear_TFI_2D
from src.msh.plot3d import PLOT3D_Block, PLOT3D
from src.msh.fluent import XF_MSH, BCType


def build_airfoil_msh(foil, ending, tk, La, Lt, R, foil_order, N):
    c0 = WingProfile(foil, ending, tk, foil_order).nurbs_rep
    c1 = Arc.from_2pnt((La, R, 0), (La, -R, 0), 180, (0, 0, 1))

    yhu = c0.start[1]
    ydu = c0.end[1]
    p = np.array([[La, yhu, 0],
                  [La, ydu, 0],
                  [La, R, 0],
                  [La, -R, 0],
                  [La + Lt, yhu, 0],
                  [La + Lt, ydu, 0],
                  [La + Lt, R, 0],
                  [La + Lt, -R, 0]])

    c2 = Line(p[0], p[2])
    c3 = Line(p[1], p[3])
    c4 = Line(p[4], p[6])
    c5 = Line(p[5], p[7])
    c6 = Line(p[4], p[0])
    c7 = Line(p[6], p[2])
    c8 = Line(p[5], p[1])
    c9 = Line(p[7], p[3])
    c10 = Line(p[1], p[0])
    c11 = Line(p[5], p[4])

    pu = [double_exponential(N[0], 0.5, -1.5, 0.5),  # c0, c1
          hyperbolic_tangent(N[1], 2),  # c2, c3, c4, c5
          single_exponential(N[2], -3),  # c6, c7, c8, c9
          np.linspace(0.0, 1.0, N[3])]  # c10, c11

    '''翼型前缘'''
    grid1 = Linear_TFI_2D(c0, c2, c1, c3)
    ppu1, ppv1 = np.meshgrid(pu[0], pu[1], indexing='ij')
    grid1.calc_grid(ppu1, ppv1)
    tm_grid1 = ThomasMiddlecoff2D(grid1.get_grid())
    tm_grid1.calc_grid()

    '''翼型后缘上部'''
    grid2 = Linear_TFI_2D(c6, c4, c7, c2)
    ppu2, ppv2 = np.meshgrid(pu[2], pu[1], indexing='ij')
    grid2.calc_grid(ppu2, ppv2)

    '''翼型后缘下部'''
    grid3 = Linear_TFI_2D(c8, c5, c9, c3)
    ppu3, ppv3 = np.meshgrid(pu[2], pu[1], indexing='ij')
    grid3.calc_grid(ppu3, ppv3)

    '''钝尾缘部分'''
    grid4 = Linear_TFI_2D(c8, c11, c6, c10)
    ppu4, ppv4 = np.meshgrid(pu[2], pu[3], indexing='ij')
    grid4.calc_grid(ppu4, ppv4)

    '''网格, 边界条件, 邻接关系'''
    blk = [tm_grid1.get_grid(), grid2.get_grid(), grid3.get_grid(), grid4.get_grid()]

    bc = [(BCType.Wall, BCType.Interior, BCType.PressureFarField, BCType.Interior),
          (BCType.Interior, BCType.PressureFarField, BCType.PressureFarField, BCType.Interior),
          (BCType.Interior, BCType.PressureFarField, BCType.PressureFarField, BCType.Interior),
          (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Wall)]

    adj = [((0, 2), (1, 4)),
           ((2, 4), (0, 4)),
           ((3, 3), (1, 1)),
           ((2, 1), (3, 1)),
           ((0, 3), (0, 0)),
           ((0, 0), (0, 1)),
           ((1, 3), (0, 0)),
           ((0, 0), (2, 3)),
           ((0, 0), (3, 4)),
           ((3, 2), (0, 0)),
           ((1, 2), (0, 0)),
           ((0, 0), (2, 2))]

    '''构建MSH文件'''
    fn = '{}_C_grid.msh'.format(foil)
    msh = XF_MSH.from_str2d_multi(blk, bc, adj)
    msh.save(fn)


class MultiBlockFluentMeshTest(unittest.TestCase):
    @staticmethod
    def test():
        build_airfoil_msh('M6', [[0, 0, 0], [10, 0, 0]], 1.2, 10, 200, 100, 5, [41, 26, 61, 2])


if __name__ == '__main__':
    unittest.main()
