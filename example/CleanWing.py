from win32com import client
import numpy as np
import math
from matplotlib import pyplot as plt
from planform import HWBWingPlanform
from load_dist import EllipticLiftDist, LinearLiftDist, UniformLiftDist
from profile import ProfileList, calc_profile_cl
from airfoil import Airfoil, airfoil_interp
from spacing import uniform, chebshev_dist_multi
from iges import IGES_Model

# catia = client.Dispatch('CATIA.Application')
# catia.Visible = True
#
# model_doc = catia.Documents
# model_part_doc = model_doc.Add('Part')
# model_part = model_part_doc.Part
# model_bodies = model_part.HybridBodies
# model_body = model_bodies.Add()
# model_body.Name = 'HWBWing'
# model_sf = model_part.HybridShapeFactory


wing_spn2 = 18.7
wing_root_len = 18.4
wing_tip_len = 1.8
wing_leading_inner_delta = (1.0, 1.55)
wing_leading_middle_delta = (1.3, 1.75)
wing_leading_outer_sweep = 28
wing_trailing_inner_delta = (0.7, -2.3)
wing_trailing_outer_spn = 13.5
wing_trailing_outer_sweep = 12
wing_planform = HWBWingPlanform(wing_root_len, wing_tip_len, wing_spn2,
                                wing_leading_inner_delta, wing_leading_middle_delta, wing_leading_outer_sweep,
                                wing_trailing_inner_delta, wing_trailing_outer_spn, wing_trailing_outer_sweep)

wing_u = uniform(51)
wing_n = len(wing_u)
wing_chord = np.array([wing_planform.chord_len(u) for u in wing_u])
wing_z = wing_u * wing_spn2

payload = 135
w = payload * 1000 * 9.8
rho = 0.4135
a = 299.5
ma = 0.8
v = ma * a
p_inf = 0.5 * rho * v ** 2

wing_lift_dist = LinearLiftDist(payload, 2 * wing_spn2, rho, v)
wing_swp_025 = np.array([wing_planform.swp_025(u) for u in wing_u])

wing_cl3 = calc_profile_cl(wing_u, wing_lift_dist, wing_planform, p_inf)
print(wing_cl3)

wing_cl2 = np.array([1.1 * wing_cl3[i] / math.cos(math.radians(wing_swp_025[i])) ** 2 for i in range(wing_n)])
print(wing_cl2)

# wing_lift_fig = plt.figure()
# wing_vc_ax = wing_lift_fig.add_subplot(111)
# fig_sp_n = 1000
# fig_x_sp = np.linspace(0, wing_spn2, fig_sp_n)
# fig_u_sp = np.linspace(0, 1, fig_sp_n)
# fig_vc_sp = [wing_lift_dist.velocity_circulation_at(u) for u in fig_u_sp]
# wing_vc_ax.plot(fig_x_sp, fig_vc_sp)
# wing_vc_ax.set_xlabel('Span-wise position')
# wing_vc_ax.set_ylabel('Velocity Circulation/({})'.format(r'$m^2 \cdot s^{-1}$'))
# wing_cl_ax = wing_vc_ax.twinx()
# wing_cl_ax.plot(wing_z, wing_cl2, '.-')
# for i in range(wing_n):
#     wing_cl_ax.text(wing_z[i], wing_cl2[i], '{:.3f}'.format(wing_cl2[i]))
# wing_cl_ax.set_ylabel('Lift coefficient')
# wing_lift_fig.tight_layout()
# wing_lift_fig.savefig('clean_wing_lift_dist.png', dpi=300)

root_airfoil = Airfoil('NACA64(3)-218A')
inner_airfoil = Airfoil('NACA63(2)-615A')
nsp = chebshev_dist_multi((0, 0.5, 1), (101, 101))
interp1, interp2 = airfoil_interp(root_airfoil, inner_airfoil, [1 / 3, 2 / 3], np.array([nsp, nsp]))
interp1.save('interp1.dat')
interp2.save('interp2.dat')

wing_foil = [root_airfoil,
             interp1,
             interp2,
             inner_airfoil,
             Airfoil('SC(2)-0714'),
             Airfoil('SC(2)-0712'),
             Airfoil('SC(2)-0712'),
             Airfoil('SC(2)-0712'),
             Airfoil('SC(2)-0712'),
             Airfoil('SC(2)-0712'),
             Airfoil('SC(2)-0712'),
             Airfoil('SC(2)-0012')]

wing_twist_ang = np.array([1.685, 0.974, 0.339, -0.358, -0.430,
                           -0.498, -0.498, -0.498, -0.498, -0.498, -0.498, 0])
wing_twist_ref = np.array([1.0] * wing_n)

wpl = ProfileList.from_planform(wing_planform, wing_foil, wing_twist_ang, wing_twist_ref, wing_u)

model = IGES_Model()
pfl = wpl.generate_nurbs_crv_list()
for c in pfl:
    model.add(c.to_iges())
model.save('HWB.igs')

# fig = plt.figure()
# ax = fig.add_subplot(111)
# wing_planform.pic(ax, u=wing_u)
# fig.tight_layout()
# fig.savefig('clean_wing.png', dpi=300)
