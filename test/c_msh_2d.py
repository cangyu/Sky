import numpy as np
import math
from src.nurbs.curve import GlobalInterpolatedCrv, Line
from src.aircraft.wing import WingProfile
from src.iges.iges_core import IGES_Model
from src.msh.spacing import single_exponential, double_exponential
from src.msh.elliptical import Laplace_2D
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
fn = "{}_C_Laplace.xyz".format(airfoil)

N0 = 35
A = 3

N1 = 40
N2 = 80
A1 = 0.5
A2 = -2.5
A3 = 0.5

La = 10
Lt = 20 * La
Lf = 15 * La
R = 10 * La
Thickness = 1.0

'''Derived Variables'''
ending = np.array([[0, 0, 0], [La, 0, 0]])
af = WingProfile(airfoil, ending, Thickness, 5)
airfoil_len = af.nurbs_rep.length()
yhu = af.pts[0][1]
ydu = af.pts[-1][1]
inner_len = 2 * Lt + airfoil_len
outer_len = 2 * Lt + ellipse_arc_len(La + Lf, R) / 2

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
c1u1 = Lt / inner_len
c1u2 = (inner_len - Lt) / inner_len


def c1(u: float):
    if u <= c1u1:
        ans = hu(u / c1u1)
    elif u <= c1u2:
        ans = af.nurbs_rep((u - c1u1) / (c1u2 - c1u1))
    else:
        ans = du((u - c1u2) / (1.0 - c1u2))

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

c3u1 = Lt / outer_len
c3u2 = (outer_len - Lt) / outer_len

'''
def c3(u: float):
    if u <= c3u1:
        ans = l0(u / c3u1)
    elif u <= c3u2:
        ans = l1((u - c3u1) / (c3u2 - c3u1))
    else:
        ans = l2((u - c3u2) / (1.0 - c3u2))

    return ans[0:2]
'''


def c3(u: float):
    return outer_crv(u)[0:2]


'''Boundary distribution'''
pv0 = single_exponential(0.0, 1.0, N0 + 1, A)
pv1 = single_exponential(0.0, 1.0, N0 + 1, A)

pu = []
pu1 = single_exponential(0.0, c1u1, N1 + 1, -A)
pu2 = double_exponential(c1u1, c1u2, N2 + 1, A1, A2, A3)
pu3 = single_exponential(c1u2, 1.0, N1 + 1, A)
for u in pu1:
    pu.append(u)

for i in range(1, len(pu2)):
    pu.append(pu2[i])

for i in range(1, len(pu3)):
    pu.append(pu3[i])

pu0 = np.copy(pu)

pu = []
pu1 = np.linspace(0.0, c3u1, N1 + 1)
pu2 = np.linspace(c3u1, c3u2, N2 + 1)
pu3 = np.linspace(c3u2, 1.0, N1 + 1)
for u in pu1:
    pu.append(u)

for i in range(1, len(pu2)):
    pu.append(pu2[i])

for i in range(1, len(pu3)):
    pu.append(pu3[i])

pu1 = np.copy(pu)

'''
grid = Laplace_2D(c1, c2, c3, c4, pu0, pu1, pv0, pv1)
grid.calc_grid()
grid.write_plot3d(fn)
'''

grid = Linear_TFI_2D(c1, c2, c3, c4)
grid.calc_msh(pu0, pu1, pv0, pv1)
grid.write_plot3d(fn)

'''Model'''
model = IGES_Model('frame.igs')
model.add_entity(af.nurbs_rep.to_iges(1, 0, [0, 0, 1]))
model.add_entity(hu.to_iges())
model.add_entity(du.to_iges())
model.add_entity(outer_crv.to_iges(1, 0, [0, 0, 1]))
model.add_entity(tul.to_iges())
model.add_entity(tdl.to_iges())
model.write()
