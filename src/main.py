from src.iges.iges_entity110 import *
from src.iges.iges_entity116 import *
from src.iges.iges_entity126 import *
from src.aircraft.wing import *
from src.nurbs.curve import *

plane = IGES_Model()

naca0012 = Airfoil("../airfoil/naca0012.dat")

z = np.array([0., 0.624029, 1.38967, 2.43503, 3.73439, 5.25574, 6.96162,
              8.81003, 10.7555, 12.75, 14.7445, 16.69, 18.5384, 20.2443, 21.7656,
              23.065, 24.1103, 24.876, 25.343, 25.5])

x_front = np.array([0., 0.05, 0.3, 1.7, 4.9, 6.85, 8.45, 9.65, 10.6, 11.1, 11.7, 12.1,
                    12.4, 12.8, 13.2, 13.7, 14.1, 14.5, 15.2, 16.])

y_front = np.linspace(0, 0, len(z))

x_tail = np.array([19.7, 19.6, 19.6, 19.5, 19.3, 19, 18.3, 17.3, 16.6, 16.5,
                   16.8, 17, 17.45, 17.8, 18.1, 18.4, 18.55, 18.65, 18.3, 17.8])

y_tail = np.linspace(0, 0, len(z))

for i in range(0, len(z)):
    plane.AddPart(IGES_Entity116(x_front[i], y_front[i], z[i]))
    plane.AddPart(IGES_Entity116(x_tail[i], y_tail[i], z[i]))

plane.AddPart(IGES_Entity110([x_front[0], y_front[0], z[0]], [x_tail[0], y_tail[0], z[0]]))
plane.AddPart(IGES_Entity110([x_front[len(z) - 1], y_front[len(z) - 1], z[len(z) - 1]], [x_tail[len(z) - 1], y_tail[len(z) - 1], z[len(z) - 1]]))

profile_pts=[]
profile_crv=[]

for i in range(0, len(z)):
    epts = np.array([[x_front[i], y_front[i], z[i]], [x_tail[i], y_tail[i], z[i]]])
    wp = Wing_Profile(naca0012, epts)
    pts = wp.getPointList()
    profile_pts.append(pts)
    '''
    for k in range(0, len(pts)):
        plane.AddPart(IGES_Entity116(pts[k][0], pts[k][1], pts[k][2]))
    '''

cc = Curve(profile_pts[0])
a, b, c = cc.generate()
ccc = IGES_Entity126(cc.p, cc.n, 1, 0, 1, 0, a, b, c, 0.0, 1.0, np.array([0,0,1.0]))
plane.AddPart(ccc)

plane.Generate()