from win32com import client
from wing import HWBWingPlanform, WingProfileList, Wing
from iges import IGES_Model
import numpy as np
from matplotlib import pyplot as plt

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

wing_foil = ['NACA0018',
             'NACA0017',
             'NACA0016',
             'NACA0015',
             'SC(2)-0614',
             'SC(2)-0612',
             'SC(2)-0712',
             'SC(2)-0712',
             'SC(2)-0712',
             'SC(2)-0712',
             'SC(2)-0612',
             'SC(2)-0412',
             'SC(2)-0012']

wing_twist_ang = np.array([0, 0, 0, 0, -2.523, -2.180, -2.637, -2.439, -2.303, -2.418, -2.306, -2.071, -2.000, ])
wing_twist_ref = np.array([1.0] * wing_n)

wpl = WingProfileList(wing_planform, wing_foil, wing_twist_ang, wing_twist_ref, wing_u)
# wing = Wing(wpl)

model = IGES_Model()
pfl = wpl.crv_list_in_nurbs()
for c in pfl:
    model.add(c.to_iges())
model.save('HWB.igs')

fig = plt.figure()
ax = fig.add_subplot(111)
wing_planform.pic(ax, u=wing_u)
fig.tight_layout()
plt.show()
