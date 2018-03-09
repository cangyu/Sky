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

'''
    Planform Design
'''
spn2 = 21.0
root_chord = 19.2
tip_chord = 1.8
leading_cpt = [(0.3, 0.8), (1.97, 2.2),
               (5.12, 4.17), (7.14, 5.47),
               (8.5, 6.95), (10.2, 9.6), (17.4, spn2)]
trailing_cpt = [(17.23, 3.83), (16.1, 6.97), (16.25, 10.15)]

_chord = [root_chord, tip_chord]
_cpt = leading_cpt + trailing_cpt
planform = HWBCommonPlanform(_chord, _cpt)

pf_pos = [0.00, 5.24, 10.48, 15.24, 20.00,
          24.54, 28.89, 34.50, 40.00,
          45.70, 55.33, 64.95, 74.57, 84.00,
          93.61, 95.87, 97.94, 100.00]
pf_u = np.array([0.01 * x for x in pf_pos])
pf_n = len(pf_u)
pf_z = np.array([x * spn2 for x in pf_u])

# Show Planform and Profiles
fig = plt.figure()
ax = fig.add_subplot(111)
planform.pic(ax, u=pf_u)
fig.tight_layout()
plt.show()

'''
    Airfoil Selection and Interpolation
'''
root_foil = Airfoil.from_local('NACA64A221')
cabin_foil = Airfoil.from_local('NACA64A221')
fusion_foil = Airfoil.from_local('NACA63A418')
inner_foil = Airfoil.from_local('NACA63A615')
wing_foil = Airfoil.from_local('SC(2)-0712')
outer_foil = Airfoil.from_local('SC(2)-0612')
tip_foil = Airfoil.from_local('SC(2)-0012')

fsp = chebshev_dist_multi([0, 0.5, 1], [101, 101])
foil_interp1 = airfoil_interp(root_foil, cabin_foil, 0.5, fsp)
foil_interp2 = airfoil_interp(cabin_foil, fusion_foil, 0.5, fsp)
foil_interp3 = airfoil_interp(fusion_foil, inner_foil, 0.5, fsp)
foil_interp4 = airfoil_interp(inner_foil, wing_foil, [1 / 3, 2 / 3], [fsp] * 2)
foil_interp5 = airfoil_interp(wing_foil, outer_foil, [1 / 5, 2 / 5, 3 / 5, 4 / 5], [fsp] * 4)
foil_interp6 = airfoil_interp(outer_foil, tip_foil, [1 / 3, 2 / 3], [fsp] * 2)

foil = [root_foil, foil_interp1, cabin_foil, foil_interp2, fusion_foil, foil_interp3, inner_foil]
foil += foil_interp4
foil.append(wing_foil)
foil += foil_interp5
foil.append(outer_foil)
foil += foil_interp6
foil.append(tip_foil)

# Save all airfoils for XFoil usage
# for k, f in enumerate(foil):
#     f.save('foil{}.dat'.format(k))

'''
    Twist and Dihedral Design
'''
param = np.array([[0.2953, 0.000, 1.00, 0.000],  # 0 - 64A221
                  [0.2998, 0.000, 1.00, 0.070],  # 1
                  [0.3006, 0.000, 1.00, 0.180],  # 2 - 64A221
                  [0.6039, 0.750, 1.00, 0.125],  # 3
                  [0.8702, 1.500, 1.00, 0.080],  # 4 - 63A418
                  [0.9114, 1.101, 1.00, 0.220],  # 5
                  [0.9855, 0.702, 1.00, 0.300],  # 6 - 63A615
                  [0.8002, 0.304, 1.00, 0.460],  # 7
                  [0.8083, -.095, 1.00, 0.575],  # 8
                  [0.7500, -.494, 1.00, 0.655],  # 9 - SC(2)-0712
                  [0.7312, -.443, 1.00, 0.705],  # 10
                  [0.7129, -.392, 1.00, 0.756],  # 11
                  [0.6939, -.341, 1.00, 0.806],  # 12
                  [0.6725, -.290, 1.00, 0.856],  # 13
                  [0.6500, -.240, 1.00, 0.905],  # 14 - SC(2)-0612
                  [0.4340, -.180, 1.00, 0.925],  # 15
                  [0.2154, -.090, 1.00, 0.940],  # 16
                  [0.0000, 0.000, 1.00, 0.945]])  # 17 - SC(2)-0012

cl2 = param[:, 0]
twist = param[:, 1]
twist_ref = param[:, 2]
y_off = param[:, 3]
wpl = ProfileList.from_planform_with_dihedral(planform, pf_u, foil, twist, twist_ref, y_off)

'''
    Write IGES file
'''
model = IGES_Model()
pfl = wpl.generate_nurbs_crv_list()
for c in pfl:
    model.add(c.to_iges())
model.save('HWB.igs')
