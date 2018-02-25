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

pf_desc = [[Airfoil.from_local('NACA64A221'), 0.00, 2],
           [Airfoil.from_local('NACA64A320'), 3.80, 2],
           [Airfoil.from_local('NACA64A420'), 8.05, 5],
           [Airfoil.from_local('NACA63A518'), 12.50, 3],
           [Airfoil.from_local('NACA63A618'), 18.24, 3],
           [Airfoil.from_local('NACA63A615'), 23.00, 12],
           [Airfoil.from_local('SC(2)-0712'), 42.64, 21],
           [Airfoil.from_local('SC(2)-0712'), 90.31, 2],
           [Airfoil.from_local('SC(2)-0612'), 93.54, 3],
           [Airfoil.from_local('SC(2)-0412'), 97.25, 3],
           [Airfoil.from_local('SC(2)-0012'), 100.00, 0]]

pf_n = len(pf_desc)
fsp = chebshev_dist_multi([0, 0.5, 1], [100, 100])

pf_pos = [pf_desc[0][1]]
foil = [pf_desc[0][0]]
for i in range(pf_n - 1):
    local_interp_pos = np.linspace(pf_desc[i][1], pf_desc[i + 1][1], pf_desc[i][2] + 2)
    pf_pos += list(local_interp_pos[1:])

    foil_interp_pos = np.linspace(0, 1, pf_desc[i][2] + 2)[1:-1]
    seg = airfoil_interp(pf_desc[i][0], pf_desc[i + 1][0], foil_interp_pos, [fsp] * len(foil_interp_pos))
    foil += seg
    foil.append(pf_desc[i + 1][0])

pf_u = np.array([0.01 * x for x in pf_pos])
pf_n = len(pf_u)
pf_z = np.array([x * spn2 for x in pf_u])

fig = plt.figure()
ax = fig.add_subplot(111)
planform.pic(ax, u=pf_u)
fig.tight_layout()
plt.show()

twist = np.zeros(pf_n)
twist_ref = np.ones(pf_n)
wpl = ProfileList.from_planform(planform, foil, twist, twist_ref, pf_u)

model = IGES_Model()
pfl = wpl.generate_nurbs_crv_list()
for c in pfl:
    model.add(c.to_iges())
model.save('HWB.igs')

# payload = 200
# rho = 0.4135  # H=10000m
# a = 299.5
# ma = 0.8
# v = ma * a
# re = 4e7
# area = planform.area
# cl = calc_global_cl(payload, area, ma, rho, a)
# print('Design Lift Coefficient: {:.4f}'.format(cl))
#
# linear_dist = LinearLiftDist(payload, planform.span, rho, v)
# elliptic_dist = EllipticLiftDist(payload, planform.span, rho, v)
# lift_dist = HybridLiftDist(payload, planform.span, rho, v)
# lift_dist.add(linear_dist, 0.35)
# lift_dist.add(elliptic_dist, 0.65)
# swp_025 = np.array([planform.swp_025(rel_pos) for rel_pos in u])
# cl3 = calc_profile_cl(u, lift_dist, planform)
# cl2 = np.array([1.1 * cl3[i] / math.cos(math.radians(swp_025[i])) ** 2 for i in range(n)])
# print(cl3)
