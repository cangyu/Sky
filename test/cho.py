import numpy as np
import math
from src.nurbs.curve import GlobalInterpolatedCrv
from src.nurbs.spline import Spline
from src.aircraft.wing import WingProfile
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from src.iges.iges_core import IGES_Model

N1 = 41
N2 = 61
La = 10
ending = np.array([[0, 0, 0], [La, 0, 0]])
Lt = 20 * La
R = 25 * La
af = WingProfile('NACA0012', ending, 1.0, 5)
hu = np.zeros((N1, 3))
du = np.zeros((N1, 3))

yhu = af.pts[0][1]
ydu = af.pts[-1][1]

for i in range(N1):
    u = i / (N1 - 1)
    hu[i][0] = (1 - u) * (La + Lt) + u * La
    hu[i][1] = yhu
    hu[i][2] = 0.0

    du[i][0] = u * (La + Lt) + (1 - u) * La
    du[i][1] = ydu
    du[i][2] = 0.0

boundary = []
for i in range(N1):
    boundary.append(hu[i])

for i in range(1, len(af.pts)):
    boundary.append(af.pts[i])

for i in range(1, N1):
    boundary.append(du[i])


def plot_line(pts):
    x = []
    y = []
    z = []
    for pnt in pts:
        x.append(pnt[0])
        y.append(pnt[1])
        z.append(pnt[2])

    mpl.rcParams['legend.fontsize'] = 10

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(x, y, z, label='parametric curve')
    ax.legend()

    plt.show()


inner = np.copy(boundary)

spl = Spline()
spl.interpolate(inner)
model = IGES_Model('inner.igs')
model.add_entity(spl.to_iges())
model.write()
