import numpy as np
import math
from matplotlib import pyplot as plt
from planform import HWBNoseBluntPlanform
from load_dist import EllipticLiftDist, LinearLiftDist, UniformLiftDist, calc_global_cl
from profile import ProfileList, calc_profile_cl, pic_profile_gamma_cl
from airfoil import Airfoil, airfoil_interp, find_alpha
from spacing import uniform, chebshev_dist_multi
from iges import IGES_Model

spn2 = 21.0
root_chord = 19.2
middle_chord = 4.50
tip_chord = 1.9
leading_cpt = [(1.2, 1.9), (2.9, 3.4), (5.9, 5.4), (8.0, 6.95), (9.3, 8.15), (10.2, 9.6), (17.4, spn2)]
trailing_cpt = [(15.15, 8.25), (15.1, 10.85)]

_chord = [root_chord, middle_chord, tip_chord]
_cpt = leading_cpt + trailing_cpt
planform = HWBNoseBluntPlanform(_chord, _cpt)

_seg = [0.] + [leading_cpt[i][1] / spn2 for i in [0, 1, 2, 5, 6]]
_num = [4, 3, 4, 5, 7]
# u = chebshev_dist_multi(_seg, _num)
u = np.array([0, 0.04, 0.09, 0.18, 0.26, 0.38, 0.51, 0.62, 0.74, 0.85, 0.93, 1])
n = len(u)
z = np.array([rel_pos * spn2 for rel_pos in u])

payload = 180
rho = 0.4135
a = 299.5
ma = 0.8
v = ma * a
re = 6e7
area = planform.area
cl_design = calc_global_cl(payload, area, ma, rho, a)
print('Design Lift Coefficient: {:.3f}'.format(cl_design))

lift_dist = LinearLiftDist(payload, planform.span, rho, v)
swp_025 = np.array([planform.swp_025(rel_pos) for rel_pos in u])
cl3 = calc_profile_cl(u, lift_dist, planform)
cl2 = np.array([1.1 * cl3[i] / math.cos(math.radians(swp_025[i])) ** 2 for i in range(n)])

# fig = plt.figure()
# vc_ax = fig.add_subplot(212)
# cl_ax = vc_ax.twinx()
# planform_ax = fig.add_subplot(211, sharex=vc_ax)
# pic_profile_gamma_cl(vc_ax, cl_ax, lift_dist, planform)
# planform.pic(planform_ax, u=u)
# fig.tight_layout()
# fig.set_size_inches(10.5, 20.5)
# fig.savefig('HWB_load_and_planform.png', dpi=300)

naca64a218 = Airfoil.from_local('NACA64(3)-218A')
naca63a615 = Airfoil.from_local('NACA63(2)-615A')
sc0714 = Airfoil.from_local('SC(2)-0714')
sc0712 = Airfoil.from_local('SC(2)-0712')
sc0012 = Airfoil.from_local('SC(2)-0012')

nsp = chebshev_dist_multi((0, 0.5, 1), (101, 101))
interp1, interp2 = airfoil_interp(naca64a218, naca63a615, [1 / 3, 2 / 3], np.array([nsp, nsp]))
# interp1.save('interp1.dat')
# interp2.save('interp2.dat')

foil = [naca64a218, interp1, interp2, naca63a615, sc0714, sc0712, sc0712, sc0712, sc0712, sc0712, sc0712, sc0012]
twist = np.array([1.685, 0.974, 0.339, -0.358, -0.430, -0.498, -0.498, -0.498, -0.498, -0.498, -0.498, 0])
# twist = np.array([find_alpha(foil[i], re, ma, cl2[i]) for i in range(n)])
twist_ref = np.ones(n)
wpl = ProfileList.from_planform(planform, foil, twist, twist_ref, u)

model = IGES_Model()
pfl = wpl.generate_nurbs_crv_list()
for c in pfl:
    model.add(c.to_iges())
model.save('HWB.igs')
