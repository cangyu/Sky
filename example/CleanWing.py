from win32com import client
import numpy as np
from matplotlib import pyplot as plt
from wing import HWBWingPlanform, WingProfileList, Wing
from wing import Airfoil, EllipticLiftDist
from grid import uniform, chebshev_dist_multi
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

wing_u = np.array([0.0000, 0.0400, 0.0985, 0.1793, 0.2600,
                   0.3850, 0.5100, 0.6200, 0.7400, 0.8500, 0.9300, 0.9800, 1.0000])
wing_n = len(wing_u)
wing_chord = np.array([wing_planform.chord_len(u) for u in wing_u])
wing_z = wing_u * wing_spn2

wing_lift_dist = EllipticLiftDist(120, wing_spn2 * 2, 0.4135, 0.75 * 299.5)

wing_cl3 = np.copy(list(map(wing_lift_dist.cl_at, wing_u, wing_chord)))
wing_cl2 = np.array([1.15 * x for x in wing_cl3])

wing_lift_fig = plt.figure()
wing_vc_ax = wing_lift_fig.add_subplot(111)
fig_sp_n = 1000
fig_x_sp = np.linspace(0, wing_spn2, fig_sp_n)
fig_u_sp = np.linspace(0, 1, fig_sp_n)
fig_vc_sp = [wing_lift_dist.velocity_circulation_at(u) for u in fig_u_sp]
wing_vc_ax.plot(fig_x_sp, fig_vc_sp)
wing_vc_ax.set_xlabel('Span-wise position')
wing_vc_ax.set_ylabel('Velocity Circulation/({})'.format(r'$m^2 \cdot s^{-1}$'))
wing_cl_ax = wing_vc_ax.twinx()
wing_cl_ax.plot(wing_z, wing_cl2, '.-')
for i in range(wing_n):
    wing_cl_ax.text(wing_z[i], wing_cl2[i], '{:.3f}'.format(wing_cl2[i]))
wing_cl_ax.set_ylabel('Lift coefficient')
wing_lift_fig.tight_layout()
wing_lift_fig.show()

wing_foil = ['NACA64(3)-218',
             'NACA0017',
             'NACA0016',
             'NACA63(2)-615_161',
             'SC(2)-0614',
             'SC(2)-0612',
             'SC(2)-0712',
             'SC(2)-0712',
             'SC(2)-0712',
             'SC(2)-0712',
             'SC(2)-0612',
             'SC(2)-0412',
             'SC(2)-0012']

wing_twist_ang = np.array([-0.202, 2, 2, -1.306,
                           -0.523, -0.180, -0.637,
                           -0.439, -0.303, -0.418, -0.306,
                           -0.071, -0.000])
wing_twist_ref = np.array([1.0] * wing_n)

wpl = WingProfileList(wing_planform, wing_foil, wing_twist_ang, wing_twist_ref, wing_u)

model = IGES_Model()
pfl = wpl.crv_list_in_nurbs()
for c in pfl:
    model.add(c.to_iges())
model.save('HWB.igs')

fig = plt.figure()
ax = fig.add_subplot(111)
wing_planform.pic(ax, u=wing_u)
fig.tight_layout()
fig.show()
