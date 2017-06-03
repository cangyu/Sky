import numpy as np
import math
from src.nurbs.curve import GlobalInterpolatedCrv, Line
from src.aircraft.wing import WingProfile
from src.iges.iges_core import IGES_Model
from src.msh.spacing import single_exponential, double_exponential
from src.msh.elliptical import Laplace_2D, Possion_2D
from src.msh.tfi import Linear_TFI_2D


def l0(u):
    return np.array([(1 - u) * (La + Lt) + u * La, R, 0])


def l1(u):
    ang = u * math.pi + math.pi / 2
    return np.array([(Lf + La) * math.cos(ang) + La, R * math.sin(ang), 0])


def l2(u):
    return np.array([(1 - u) * La + u * (La + Lt), -R, 0])


def ellipse_arc_len(a, b):
    if a < b:
        a, b = b, a

    return math.pi * (3 * (a + b) - math.sqrt((3 * a + b) * (a + 3 * b)))


'''Design Variables'''
airfoil = 'NACA0012'

N0 = 25
A = 2

N1 = 30
N2 = 70
A1 = 0.5
A2 = -2.5
A3 = 0.5
A4 = 2.5

La = 10
Lt = 20 * La
Lf = 15 * La
R = 10 * La
Thickness = 1.0

tfi = True
laplace = True
possion = False
build_model = False

'''Derived Variables'''
ending = np.array([[0, 0, 0], [La, 0, 0]])
af = WingProfile(airfoil, ending, Thickness, 5)
yhu = af.pts[0][1]
ydu = af.pts[-1][1]
T1 = 2 * Lt + af.nurbs_rep.length()
T3 = 2 * Lt + ellipse_arc_len(La + Lf, R) / 2
u1 = 1 / 3
u2 = 2 / 3
# u1 = Lt / T1
# u2 = (T1 - Lt) / T1

# u1 = Lt / T3
# u2 = (T3 - Lt) / T3

'''C2,C4'''
tul = Line([La + Lt, yhu, 0], [La + Lt, R, 0])
tdl = Line([La + Lt, ydu, 0], [La + Lt, -R, 0])


def c2(u):
    return tul(u)[0:2]


def c4(u):
    return tdl(u)[0:2]


'''C1'''
hu = Line([La + Lt, yhu, 0], [La, yhu, 0])
du = Line([La, ydu, 0], [La + Lt, ydu, 0])


def c1(u: float):
    if u <= u1:
        ans = hu(u / u1)
    elif u <= u2:
        ans = af.nurbs_rep((u - u1) / (u2 - u1))
    else:
        ans = du((u - u2) / (1.0 - u2))

    return ans[0:2]


'''C3'''
outer = []
for i in range(N1):
    outer.append(l0(i / (N1 - 1)))

for i in range(1, N2):
    outer.append(l1(i / (N2 - 1)))

for i in range(1, N1):
    outer.append(l2(i / (N1 - 1)))

outer = np.copy(outer)
outer_crv = GlobalInterpolatedCrv(outer, method='chord')


def c3(u: float):
    if u <= u1:
        ans = l0(u / u1)
    elif u <= u2:
        ans = l1((u - u1) / (u2 - u1))
    else:
        ans = l2((u - u2) / (1.0 - u2))

    return ans[0:2]


'''
def c3(u: float):
    return outer_crv(u)[0:2]
'''

'''Boundary distribution'''
pv = single_exponential(0.0, 1.0, N0 + 1, A)
# pv = np.linspace(0.0, 1.0, N0 + 1)
# pu = np.linspace(0.0, 1.0, 2 * N1 + N2 + 1)
# pu = double_exponential(0.0, 1.0, 2 * N1 + N2 + 1, A1, A2, A3)


pu = []
pu1 = single_exponential(0.0, u1, N1 + 1, -A4)
pu2 = double_exponential(u1, u2, N2 + 1, A1, A2, A3)
pu3 = single_exponential(u2, 1.0, N1 + 1, A4)
for u in pu1:
    pu.append(u)

for i in range(1, len(pu2)):
    pu.append(pu2[i])

for i in range(1, len(pu3)):
    pu.append(pu3[i])

pu = np.copy(pu)

if tfi:
    grid = Linear_TFI_2D(c1, c2, c3, c4)
    grid.calc_msh(pu, pv)
    grid.write_plot3d(airfoil + "_TFI.xyz")

if laplace:
    grid = Laplace_2D(c1, c2, c3, c4, pu, pv)
    grid.calc_grid()
    grid.write_plot3d(airfoil + "_Laplace.xyz")

if possion:
    grid = Possion_2D(c1, c2, c3, c4, pu, pv)
    grid.calc_grid(100)
    grid.write_plot3d(airfoil + "_Possion.xyz")

if build_model:
    model = IGES_Model('frame.igs')
    model.add_entity(af.nurbs_rep.to_iges)
    model.add_entity(hu.to_iges)
    model.add_entity(du.to_iges)
    model.add_entity(outer_crv.to_iges(1, 0, [0, 0, 1]))
    model.add_entity(tul.to_iges)
    model.add_entity(tdl.to_iges)
    model.write()
