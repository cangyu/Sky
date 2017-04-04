import numpy as np
from src.iges.iges_entity110 import *
from src.iges.iges_entity114 import *
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

'''
plane.AddPart(IGES_Entity112_Builder(z, x_front, y_front, z).GetEntity())
plane.AddPart(IGES_Entity112_Builder(z, x_tail, y_tail, z).GetEntity())

plane.AddPart(IGES_Entity110_Builder(root_line).GetEntity())
plane.AddPart(IGES_Entity110_Builder(tip_line).GetEntity())

for i in range(0, len(z)):
    plane.AddPart(IGES_Entity116_Builder(x_front[i], y_front[i], z[i]).GetEntity())
    plane.AddPart(IGES_Entity116_Builder(x_tail[i], y_tail[i], z[i]).GetEntity())
'''

naca0012 = Airfoil("../airfoil/naca0012.dat")

profile = []

for i in range(0, 1):
    epts = np.array([[x_front[i], y_front[i], z[i]],
                     [x_tail[i], y_tail[i], z[i]]])

    wp = Wing_Profile(naca0012, epts)
    profile.append(wp)
    wp.AttachTo(plane)

plane.Generate()


u = np.zeros(len(z), dtype=float)
for i in range(0, len(z)):
    u[i] = z[i]

v = np.zeros(len(naca0012.x))
for i in range(0, len(naca0012.x)):
    v[i] = naca0012.x[i]

xx = np.zeros((len(u), len(v)))
yy = np.zeros((len(u), len(v)))
zz = np.zeros((len(u), len(v)))

for i in range(0, len(u)):
    for j in range(0, len(v)):
        xx[i][j] = profile[i].pts[0][0][j]
        yy[i][j] = profile[i].pts[0][1][j]
        zz[i][j] = profile[i].pts[0][2][j]

fx = interpolate.RectBivariateSpline(u, v, xx)
fy = interpolate.RectBivariateSpline(u, v, yy)
fz = interpolate.RectBivariateSpline(u, v, zz)
f = [fx, fy, fz]



