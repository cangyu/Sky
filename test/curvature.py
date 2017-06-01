import numpy as np
from src.aircraft.wing import Airfoil
from src.nurbs.curve import GlobalInterpolatedCrv
import matplotlib.pyplot as plt

L = 100
af = Airfoil('M6')
for i in range(len(af.pts)):
    af.pts[i] *= L

crv = GlobalInterpolatedCrv(af.pts, 5, 'centripetal')

N = 2000
u = np.linspace(0, 1.0, N)
kappa = np.zeros(N)
for i in range(N):
    kappa[i] = crv.curvature(u[i])

print(kappa)

plt.figure()
plt.plot(u, kappa)
plt.show()
