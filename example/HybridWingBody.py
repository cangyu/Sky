import numpy as np
import math
from matplotlib import pyplot as plt
from planform import HWBCommonPlanform
from airfoil import Airfoil, airfoil_interp
from profile import ProfileList, calc_profile_cl
from iges import IGES_Model
from load_dist import calc_global_cl
from load_dist import LinearLiftDist, EllipticLiftDist, HybridLiftDist
from spacing import chebshev_dist_multi

spn2 = 21.0
root_chord = 19.2
middle_chord = 4.50
tip_chord = 2.0
leading_cpt = [(0.18, 0.8), (1.97, 2.2),
               (5.12, 4.17), (7.14, 5.47),
               (8.58, 6.95), (10.2, 9.6), (17.4, spn2)]
trailing_cpt = [(17.23, 3.83), (16.1, 6.97), (16.3, 10.15)]

_chord = [root_chord, tip_chord]
_cpt = leading_cpt + trailing_cpt
planform = HWBCommonPlanform(_chord, _cpt)

pos = [0, 1.91, 3.81, 5.92, 8.03,
       10.48, 14.52, 18.24, 22.86, 26.05,
       29.85, 33.10, 36.91, 40.72, 44.56,
       48.33, 54.79, 61.25, 67.71, 74.17,
       80.62, 87.08, 90.31, 93.54, 97.25, 100]
u = np.array([0.01 * x for x in pos])
n = len(u)
z = np.array([x * spn2 for x in u])

fig = plt.figure()
ax = fig.add_subplot(111)
planform.pic(ax, u=u)
fig.tight_layout()
plt.show()

payload = 200
rho = 0.4135  # H=10000m
a = 299.5
ma = 0.8
v = ma * a
re = 4e7
area = planform.area
cl = calc_global_cl(payload, area, ma, rho, a)
print('Design Lift Coefficient: {:.4f}'.format(cl))

linear_dist = LinearLiftDist(payload, planform.span, rho, v)
elliptic_dist = EllipticLiftDist(payload, planform.span, rho, v)
lift_dist = HybridLiftDist(payload, planform.span, rho, v)
lift_dist.add(linear_dist, 0.35)
lift_dist.add(elliptic_dist, 0.65)
swp_025 = np.array([planform.swp_025(rel_pos) for rel_pos in u])
cl3 = calc_profile_cl(u, lift_dist, planform)
cl2 = np.array([1.1 * cl3[i] / math.cos(math.radians(swp_025[i])) ** 2 for i in range(n)])
print(cl3)

fsp = chebshev_dist_multi([0, 0.5, 1], [100, 100])

naca64a221 = Airfoil.from_local('NACA64A221')
naca64a420 = Airfoil.from_local('NACA64A420')
foil_interp1, foil_interp2 = airfoil_interp(naca64a221, naca64a420, [1 / 3, 2 / 3], [fsp, fsp])
naca64a618 = Airfoil.from_local('NACA64A618')
foil_interp3, foil_interp4 = airfoil_interp(naca64a420, naca64a618, [1 / 3, 2 / 3], [fsp, fsp])
sc0712 = Airfoil.from_local('SC(2)-0712')
middle_interp = airfoil_interp(naca64a618, sc0712, [1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6], [fsp] * 5)
sc0412 = Airfoil.from_local('SC(2)-0412')
foil_interp5 = airfoil_interp(sc0712, sc0412, 0.5, fsp)
sc0012 = Airfoil.from_local('SC(2)-0012')
foil_interp6 = airfoil_interp(sc0412, sc0012, 0.5, fsp)

foil = [naca64a221, naca64a221, naca64a221,
        foil_interp1, foil_interp2, naca64a420,
        foil_interp3, foil_interp4, naca64a618,
        middle_interp[0], middle_interp[1], middle_interp[2], middle_interp[3], middle_interp[4],
        sc0712, sc0712, sc0712, sc0712, sc0712, sc0712, sc0712, sc0712,
        foil_interp5, sc0412, foil_interp6, sc0012]
twist = np.zeros(n)
twist_ref = np.ones(n)
wpl = ProfileList.from_planform(planform, foil, twist, twist_ref, u)

model = IGES_Model()
pfl = wpl.generate_nurbs_crv_list()
for c in pfl:
    model.add(c.to_iges())
model.save('HWB.igs')
