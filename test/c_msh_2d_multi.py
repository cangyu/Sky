import numpy as np
import math
from src.nurbs.curve import GlobalInterpolatedCrv, Line
from src.aircraft.wing import WingProfile
from src.msh.spacing import single_exponential, double_exponential
from src.msh.elliptical import Laplace_2D, Possion_2D
from src.msh.tfi import Linear_TFI_2D
from src.msh.plot3d import Plot3D_MultiBlock

airfoil = 'M6'
Thickness = 1.0
La = 10
Lt = 20 * La
Lf = 15 * La
R = 10 * La
N = np.array([81, 46, 91, 2])

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
    return C0(u)[0:2]


def c1(u):
    return C1(u)[0:2]


def c2(u):
    return C2(u)[0:2]


def c3(u):
    return C3(u)[0:2]


def c4(u):
    return C4(u)[0:2]


def c5(u):
    return C5(u)[0:2]


def c6(u):
    return C6(u)[0:2]


def c7(u):
    return C7(u)[0:2]


def c8(u):
    return C8(u)[0:2]


def c9(u):
    return C9(u)[0:2]


def c10(u):
    return C10(u)[0:2]


def c11(u):
    return C11(u)[0:2]


pu = []
pu.append(double_exponential(0.0, 1.0, N[0], 0.5, -1.5, 0.5))
pu.append(single_exponential(0.0, 1.0, N[1], 3))
pu.append(single_exponential(0.0, 1.0, N[2], -3))
pu.append(np.linspace(0.0, 1.0, N[3]))

blk_list = []
grid = Laplace_2D(c0, c2, c1, c3, pu[0], pu[1], 1.0, 1.0)
grid.calc_grid()
blk_list.append(grid.plot3d_blk())

grid = Linear_TFI_2D(c8, c11, c6, c10)
grid.calc_msh(pu[2], pu[3])
blk_list.append(grid.plot3d_blk())

grid = Linear_TFI_2D(c7, c4, c7, c2)
grid.calc_msh(pu[2], pu[1])
blk_list.append(grid.plot3d_blk())

grid = Linear_TFI_2D(c8, c5, c9, c3)
grid.calc_msh(pu[2], pu[1])
blk_list.append(grid.plot3d_blk())

msh = Plot3D_MultiBlock(blk_list)
msh.output('{}_C_grid.xyz'.format(airfoil))
