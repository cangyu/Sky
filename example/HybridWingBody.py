import numpy as np
import math
from matplotlib import pyplot as plt
from planform import HWBCommonPlanform
from airfoil import Airfoil, airfoil_interp
from profile import ProfileList
from iges import IGES_Model
from spacing import chebshev_dist_multi

'''
    Planform Design
'''
fuselage_len = 27.0
spn2 = 21.0
root_chord = 19.2
tip_chord = 2.0
leading_cpt = [(0.3, 0.8), (1.97, 2.2),
               (5.12, 4.17), (7.00, 5.47),
               (8.45, 6.95), (9.91, 9.6), (15.20, spn2)]
trailing_cpt = [(17.1, 3.83), (16.0, 6.97), (15.95, 10.15)]

_chord = [root_chord, tip_chord]
_cpt = leading_cpt + trailing_cpt
planform = HWBCommonPlanform(_chord, _cpt)
mac = planform.mean_aerodynamic_chord
area2 = planform.area / 2
print('Planform Info:')
print('  MAC: {:.3f} m'.format(mac))
print('  Area: {:.2f} m^2'.format(area2))

pf_pos = [0.00, 10.48, 20.00, 26.55, 33.10, 39.40,
          45.70, 61.35, 77.20, 93.50, 96.85, 100.00]
pf_u = np.array([0.01 * x for x in pf_pos])
pf_n = len(pf_u)
pf_z = np.array([x * spn2 for x in pf_u])
pf_swp = np.array([planform.swp_025(_u) for _u in pf_u])
print('25% Sweep on each profile:')
for i in range(pf_n):
    print('{:>3d}: {:.2f}'.format(i, pf_swp[i]))

# Show Planform and Profiles
fig = plt.figure()
ax = fig.add_subplot(111)
planform.pic(ax, u=pf_u)
lx = np.array([x[1] for x in leading_cpt])
ly = np.array([x[0] for x in leading_cpt])
ax.scatter(lx, ly, c='r', marker='x')
tx = np.array([x[1] for x in trailing_cpt])
ty = np.array([x[0] for x in trailing_cpt])
ax.scatter(tx, ty, c='r', marker='x')
fig.tight_layout()
plt.show()

'''
    Airfoil Selection and Interpolation
'''
root_foil = Airfoil.from_local('NACA64A021')
cabin_foil = Airfoil.from_local('NACA64A220')
fusion_foil = Airfoil.from_local('SC(2)-0518')
inner_foil = Airfoil.from_local('SC(2)-0614')
wing_foil = Airfoil.from_local('SC(2)-0712')
outer_foil = Airfoil.from_local('SC(2)-0612')
extend_foil = Airfoil.from_local('SC(2)-0412')
tip_foil = Airfoil.from_local('SC(2)-0012')

fsp = chebshev_dist_multi([0, 0.5, 1], [101, 101])
foil_interp3 = airfoil_interp(fusion_foil, inner_foil, 0.5, fsp)
foil_interp4 = airfoil_interp(inner_foil, wing_foil, 0.5, fsp)
foil_interp5 = airfoil_interp(wing_foil, outer_foil, [1 / 3, 2 / 3], [fsp] * 2)

foil = [root_foil, cabin_foil, fusion_foil, foil_interp3, inner_foil, foil_interp4, wing_foil]
foil += foil_interp5
foil.append(outer_foil)
foil.append(extend_foil)
foil.append(tip_foil)

# Save all airfoils for XFoil usage
for k, f in enumerate(foil):
    f.save('foil{}.dat'.format(k))

'''
    Twist and Dihedral Design
'''
Ma = 0.8
a = 299.5
kinematic_viscosity = 3.52509e-5
Re = a * Ma * mac / kinematic_viscosity
print('Mach Num: {:.2f}'.format(Ma))
print('Reynolds Num: {:.2f}'.format(Re))

pf_ma = [Ma * math.cos(math.radians(pf_swp[i])) for i in range(pf_n)]
print('Ma_2D on each profile:')
for i in range(pf_n):
    print('  {:>2d}: {:.4f}'.format(i, pf_ma[i]))

pf_tc = np.array([0.21, 0.2, 0.18, 0.16, 0.14, 0.13,
                  0.12, 0.12, 0.12, 0.12, 0.12, 0.12])
print('Relative thickness on each profile:')
for i in range(pf_n):
    print('  {:>2d}: {:.2f}'.format(i, pf_tc[i] * 100))

estimated_cl = [10 * (0.95 - pf_tc[i] - pf_ma[i]) for i in range(pf_n)]
print('Estimated CL_2D on each profile:')
for i in range(pf_n):
    print('  {:>2d}: {:.4f}'.format(i, estimated_cl[i]))

# Profile Params:   Cl_2D  Twist    Pos  Y_off
param = np.array([[0.0000, 0.500, 0.375, 0.000],  # 0 - 64A021
                  [0.4148, 2.000, 0.380, 0.145],  # 1 - 64A220
                  [0.6409, 1.500, 0.480, 0.730],  # 2 - SC(2)-0518
                  [0.0000, 0.850, 0.480, 1.050],  # 3
                  [0.6199, -.150, 1.000, 1.240],  # 4 - SC(2)-0614
                  [0.0000, -.285, 1.000, 1.380],  # 5
                  [0.7018, -.435, 1.000, 1.490],  # 6 - SC(2)-0712
                  [0.0000, -.350, 1.000, 1.770],  # 7
                  [0.0000, -.270, 1.000, 2.020],  # 8
                  [0.6057, -.185, 1.000, 2.220],  # 9 - SC(2)-0612
                  [0.3400, -.025, 1.000, 2.300],  # 10 - SC(2)-0412
                  [0.0000, 0.000, 1.000, 2.335]])  # 11 - SC(2)-0012

wpl = ProfileList.from_planform_with_dihedral(planform, pf_u, foil, param[:, 1], param[:, 2], param[:, 3])

'''
    Write IGES file
'''
model = IGES_Model()
pfl = wpl.generate_nurbs_crv_list()
for c in pfl:
    model.add(c.to_iges())
model.save('HWB.igs')
