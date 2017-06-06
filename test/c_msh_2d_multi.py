import numpy as np
import math
from src.nurbs.curve import GlobalInterpolatedCrv, Line
from src.aircraft.wing import WingProfile
from src.msh.spacing import single_exponential, double_exponential, hyperbolic_tangent
from src.msh.elliptical import Laplace2D, ThomasMiddlecoff2D
from src.msh.tfi import Linear_TFI_2D
from src.msh.plot3d import PLOT3D_Block, PLOT3D

airfoil = 'M6'
Thickness = 1.0
La = 10
Lt = 20 * La
Lf = 15 * La
R = 10 * La
N = np.array([41, 26, 61, 2])

ending = np.array([[0, 0, 0], [La, 0, 0]])
inner = WingProfile(airfoil, ending, Thickness, 5)
yhu = inner.pts[0][1]
ydu = inner.pts[-1][1]

p = np.array([[La, yhu, 0],
              [La, ydu, 0],
              [La, R, 0],
              [La, -R, 0],
              [La + Lt, yhu, 0],
              [La + Lt, ydu, 0],
              [La + Lt, R, 0],
              [La + Lt, -R, 0]])


def ellipse_arc_len(a, b):
    if a < b:
        a, b = b, a

    return math.pi * (3 * (a + b) - math.sqrt((3 * a + b) * (a + 3 * b)))


def ang(u):
    return u * math.pi + math.pi / 2


def l1(u):
    return np.array([(Lf + La) * math.cos(ang(u)) + La, R * math.sin(ang(u)), 0])


outer = np.zeros((N[0], 3))
for i in range(N[0]):
    outer[i] = np.copy(l1(i / (N[0] - 1)))

C0 = inner.nurbs_rep
C1 = GlobalInterpolatedCrv(outer, method='chord')
C2 = Line(p[0], p[2])
C3 = Line(p[1], p[3])
C4 = Line(p[4], p[6])
C5 = Line(p[5], p[7])
C6 = Line(p[4], p[0])
C7 = Line(p[6], p[2])
C8 = Line(p[5], p[1])
C9 = Line(p[7], p[3])
C10 = Line(p[1], p[0])
C11 = Line(p[5], p[4])


def c0(u):
    return C0(u)


def c1(u):
    return C1(u)


def c2(u):
    return C2(u)


def c3(u):
    return C3(u)


def c4(u):
    return C4(u)


def c5(u):
    return C5(u)


def c6(u):
    return C6(u)


def c7(u):
    return C7(u)


def c8(u):
    return C8(u)


def c9(u):
    return C9(u)


def c10(u):
    return C10(u)


def c11(u):
    return C11(u)


pu = [double_exponential(N[0], 0.5, -1.5, 0.5),  # c0, c1
      hyperbolic_tangent(N[1], 2),  # c2, c3, c4, c5
      single_exponential(N[2], -3),  # c6, c7, c8, c9
      np.linspace(0.0, 1.0, N[3])]  # c10, c11

msh = PLOT3D()
"""翼型前缘"""
grid1 = Linear_TFI_2D(c0, c2, c1, c3)
ppu1, ppv1 = np.meshgrid(pu[0], pu[1], indexing='ij')
grid1.calc_grid(ppu1, ppv1)
lgrid1 = Laplace2D(grid1.get_grid())
lgrid1.calc_grid()
g1blk = PLOT3D_Block.build_from_2d(lgrid1.get_grid())
g1blk.set_area_iblank([0, N[0] - 1], range(1, N[1] - 1), [0], -2)  # c2,c3
msh.add_block(g1blk)

"""翼型后缘上部"""
grid2 = Linear_TFI_2D(c6, c4, c7, c2)
ppu2, ppv2 = np.meshgrid(pu[2], pu[1], indexing='ij')
grid2.calc_grid(ppu2, ppv2)
g2blk = PLOT3D_Block.build_from_2d(grid2.get_grid())
g2blk.set_area_iblank([N[2] - 1], range(1, N[1] - 1), [0], -1)  # c2
g2blk.set_area_iblank(range(1, N[2] - 1), [0], [0], -4)  # c6
msh.add_block(g2blk)

"""翼型后缘下部"""
grid3 = Linear_TFI_2D(c8, c5, c9, c3)
ppu3, ppv3 = np.meshgrid(pu[2], pu[1], indexing='ij')
grid3.calc_grid(ppu3, ppv3)
g3blk = PLOT3D_Block.build_from_2d(grid3.get_grid())
g2blk.set_area_iblank([N[2] - 1], range(1, N[1] - 1), [0], -1)  # c3
g2blk.set_area_iblank(range(1, N[2] - 1), [0], [0], -4)  # c8
msh.add_block(g3blk)

"""钝尾缘部分"""
grid4 = Linear_TFI_2D(c8, c11, c6, c10)
ppu4, ppv4 = np.meshgrid(pu[2], pu[3], indexing='ij')
grid4.calc_grid(ppu4, ppv4)
g4blk = PLOT3D_Block.build_from_2d(grid4.get_grid())
g4blk.set_area_iblank(range(1, N[2] - 1), [0], [0], -3)  # c8
g4blk.set_area_iblank(range(1, N[2] - 1), [N[3] - 1], [0], -2)  # c6
msh.add_block(g4blk)
msh.write('{}_C_grid.xyz'.format(airfoil))
