import numpy as np
from src.iges.iges_entity110 import *
from src.iges.iges_entity116 import *
from src.wing import *

z = np.array([0., 0.624029, 1.38967, 2.43503, 3.73439, 5.25574, 6.96162,
              8.81003, 10.7555, 12.75, 14.7445, 16.69, 18.5384, 20.2443, 21.7656,
              23.065, 24.1103, 24.876, 25.343, 25.5])

x_front = np.array([0., 0.05, 0.3, 1.7, 4.9, 6.85, 8.45, 9.65, 10.6, 11.1, 11.7, 12.1,
                    12.4, 12.8, 13.2, 13.7, 14.1, 14.5, 15.2, 16.])

y_front = np.linspace(0, 0, len(z))

x_tail = np.array([19.7, 19.6, 19.6, 19.5, 19.3, 19, 18.3, 17.3, 16.6, 16.5,
                   16.8, 17, 17.45, 17.8, 18.1, 18.4, 18.55, 18.65, 18.3, 17.8])

y_tail = np.linspace(0, 0, len(z))

root_line = np.array([[x_front[0], y_front[0], z[0]],
                      [x_tail[0], y_tail[0], z[0]]])

tip_line = np.array([[x_front[len(z) - 1], y_front[len(z) - 1], z[len(z) - 1]],
                     [x_tail[len(z) - 1], y_tail[len(z) - 1], z[len(z) - 1]]])

plane = IGES_Model()

plane.AddPart(IGES_Entity112_Builder(z, x_front, y_front, z).GetEntity())
plane.AddPart(IGES_Entity112_Builder(z, x_tail, y_tail, z).GetEntity())

plane.AddPart(IGES_Entity110_Builder(root_line).GetEntity())
plane.AddPart(IGES_Entity110_Builder(tip_line).GetEntity())

for i in range(0, len(z)):
    plane.AddPart(IGES_Entity116_Builder(x_front[i], y_front[i], z[i]).GetEntity())
    plane.AddPart(IGES_Entity116_Builder(x_tail[i], y_tail[i], z[i]).GetEntity())

naca0012 = Airfoil("../airfoil/naca0012.dat")

profile = []

for i in range(0, len(z)):
    epts = np.array([[x_front[i], y_front[i], z[i]],
                     [x_tail[i], y_tail[i], z[i]]])

    wp = Wing_Profile(naca0012, epts)
    profile.append(wp)
    wp.AttachTo(plane)

plane.Generate()

from scipy import interpolate

u = np.zeros(len(z), dtype=float)
for i in range(0, len(z)):
    u[i] = z[i]

v = np.zeros(len(naca0012.x))
for i in range(0, len(naca0012.x)):
    v[i] = naca0012.x[i]

xx = np.zeros((len(v), len(u)))
yy = np.zeros((len(v), len(u)))
zz = np.zeros((len(v), len(u)))

for i in range(0, len(v)):
    for j in range(0, len(u)):
        xx[i][j] = profile[j].pts[0][0][i]
        yy[i][j] = profile[j].pts[0][1][i]
        zz[i][j] = profile[j].pts[0][2][i]

fx = interpolate.interp2d(u, v, xx, kind='cubic')
fy = interpolate.interp2d(u, v, yy, kind='cubic')
fz = interpolate.interp2d(u, v, zz, kind='cubic')

uu = np.arange(0, 25.5, 0.1)
vv = np.arange(0, 1, 0.05)

nx = fx(uu, vv)
ny = fy(uu, vv)
nz = fz(uu, vv)

import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

fig = plt.figure()
ax = fig.gca(projection='3d')

uuu, vvv = np.meshgrid(uu, vv)

# Plot the surface.
surf = ax.plot_surface(uuu, vvv, ny, cmap=cm.coolwarm, linewidth=0, antialiased=True)
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
