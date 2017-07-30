import numpy as np
import math
from copy import deepcopy
from src.aircraft.wing import Wing, WingProfile
from src.iges.iges_core import IGES_Model
from src.nurbs.curve import Arc, Line, GlobalInterpolatedCrv
from src.nurbs.surface import ExtrudedSurf, RuledSurf, Coons

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

fn = "GroundEffectCraft.igs"
gec = IGES_Model(fn)

nose_len = 3.3
fuselage_diameter = 2.95
fuselage_len = 12.5
tail_end_diameter = 0.4
tail_tip_back_angle = 25
tail_y_down_delta = 0.3
wing_x_offset = 3.2
wing_y_offset = -0.4
vs_x_offset = 0.8
vs_y_offset = 1.4
hs_x_offset = 2.3
hs_y_offset = 5.1

'''Nose'''
nose_x = np.array([0, 0.4, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
nose_up_y = np.array([0.05, 1.6, 2.25, 3, 3.3, 3.55, 4, 5, 5.5, 5.82, 6, 6.1])
nose_bottom_y = np.array([-0.1, -1.6, -2, -2.55, -2.8, -3, -3, -3.05, -3.05, -3.1, -3.1, -3.2])
nose_side_z = np.array([0.2, 2, 2.8, 3.4, 4, 4.2, 4.4, 4.48, 4.5, 4.55, 4.59, 4.6])
nose_side_y = np.copy(nose_x)

n = len(nose_x)
nose_height = nose_up_y[-1] - nose_bottom_y[-1]
fuselage_radius = fuselage_diameter / 2
fuselage_center_y = 0.5 * (nose_up_y[-1] + nose_bottom_y[-1]) / nose_height * fuselage_diameter
for k in range(n):
    nose_x[k] /= nose_x[-1]
    nose_x[k] *= nose_len
    nose_up_y[k] /= nose_height
    nose_up_y[k] *= fuselage_diameter
    nose_bottom_y[k] /= nose_height
    nose_bottom_y[k] *= fuselage_diameter
    nose_side_z[k] /= nose_side_z[-1]
    nose_side_z[k] *= fuselage_radius
    nose_side_y[k] /= nose_side_y[-1]
    nose_side_y[k] *= fuselage_center_y

nose_pts = np.empty((3, n, 3), float)
for k in range(n):
    nose_pts[0][k] = np.array([nose_x[k], nose_up_y[k], 0])
    nose_pts[1][k] = np.array([nose_x[k], nose_bottom_y[k], 0])
    nose_pts[2][k] = np.array([nose_x[k], nose_side_y[k], nose_side_z[k]])

nose_c1 = GlobalInterpolatedCrv(nose_pts[0])
nose_c2 = GlobalInterpolatedCrv(nose_pts[1])
nose_c3 = GlobalInterpolatedCrv(nose_pts[2])
nose_arc1 = Arc.from_2pnt(nose_pts[0][0], nose_pts[2][0], 90, (1, 0, 0))
nose_arc2 = Arc.from_2pnt(nose_pts[2][0], nose_pts[1][0], 90, (1, 0, 0))
nose_arc3 = Arc.from_2pnt(nose_pts[0][-1], nose_pts[2][-1], 90, (1, 0, 0))
nose_arc4 = Arc.from_2pnt(nose_pts[2][-1], nose_pts[1][-1], 90, (1, 0, 0))

nose_surf1 = Coons(nose_arc1, nose_arc3, nose_c1, nose_c3)
nose_surf2 = Coons(nose_arc2, nose_arc4, nose_c3, nose_c2)
ns1m = deepcopy(nose_surf1)
ns2m = deepcopy(nose_surf2)
ns1m.mirror('Z')
ns2m.mirror('Z')

gec.add_entity(nose_surf1.to_iges())
gec.add_entity(nose_surf2.to_iges())
gec.add_entity(ns1m.to_iges())
gec.add_entity(ns2m.to_iges())

'''Fuselage'''
fuselage_x_start = nose_len
body_profile = Arc.from_2pnt(nose_pts[0][-1], nose_pts[1][-1], 180, (1, 0, 0))
body_half = ExtrudedSurf(body_profile, (fuselage_len, 0, 0))
bhm = deepcopy(body_half)
bhm.mirror('Z')

gec.add_entity(body_half.to_iges())
gec.add_entity(bhm.to_iges())

'''Tail'''
tail_len = (fuselage_diameter - tail_end_diameter) / math.tan(math.radians(tail_tip_back_angle))
tail_x_start = fuselage_x_start + fuselage_len
tail_x_end = tail_x_start + tail_len
p1 = np.array([tail_x_start, fuselage_center_y + fuselage_radius, 0])
p2 = np.array([tail_x_start, fuselage_center_y - fuselage_radius, 0])
p3 = np.array([tail_x_end, p1[1] - tail_y_down_delta, 0])
p4 = np.array([tail_x_end, p3[1] - tail_end_diameter, 0])
tail_circle1 = Arc.from_2pnt(p1, p2, 180, (1, 0, 0))
tail_circle2 = Arc.from_2pnt(p3, p4, 180, (1, 0, 0))
tail_top_line = Line(p1, p3)
tail_bottom_line = Line(p2, p4)
tail_surf = Coons(tail_circle1, tail_circle2, tail_top_line, tail_bottom_line)
tsm = deepcopy(tail_surf)
tsm.mirror('Z')

gec.add_entity(tail_surf.to_iges())
gec.add_entity(tsm.to_iges())

'''Wing'''
wing_x_start = fuselage_x_start + wing_x_offset
wing_pan_dir = (wing_x_start, wing_y_offset, 0)

foil = ['NACA0012', 'M6', 'M6']
z_offset = np.array([0, 7.5, 9.2])
length = np.array([3.3, 3.3, 1.6])
sweep_back = np.array([0, 0, 3], float)
twist = np.array([0, 0, 3], float)
dihedral = np.array([0, 3, 5], float)
twist_pos = np.array([0.25, 0.25, 0.25])
y_ref = np.zeros(3)
thickness_factor = np.array([1.1, 1.0, 0.9], float)

crv_list = []
for k in range(3):
    wp = WingProfile.from_geom_param(foil[k], z_offset[k], length[k], sweep_back[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness_factor[k])
    crv = wp.nurbs_rep()
    crv.pan(wing_pan_dir)
    crv_list.append(crv)

inner_surf = RuledSurf(crv_list[0], crv_list[1])
outer_surf = RuledSurf(crv_list[1], crv_list[2])

gec.add_entity(inner_surf.to_iges())
gec.add_entity(outer_surf.to_iges())

ism = deepcopy(inner_surf)
osm = deepcopy(outer_surf)
ism.mirror('Z')
osm.mirror('Z')

gec.add_entity(ism.to_iges())
gec.add_entity(osm.to_iges())

'''Vertical Stablizer'''
vs_x_start = tail_x_start + vs_x_offset
vs_y_start = vs_y_offset
vs_pan_dir = (vs_x_start, vs_y_start, 0)

foil = ['NACA0012', 'NACA0012']
z_offset = np.array([0, 3.7])
length = np.array([3.1, 1.5])
sweep_back = np.array([0, 25.0])
twist = np.zeros(2)
dihedral = np.zeros(2)
twist_pos = np.array([0.25, 0.25])
y_ref = np.zeros(2)
thickness_factor = np.ones(2) * 2

crv_list = []
for k in range(2):
    wp = WingProfile.from_geom_param(foil[k], z_offset[k], length[k], sweep_back[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness_factor[k])
    crv = wp.nurbs_rep()
    crv.rotate((0, 0, 0), (-1, 0, 0), 90)
    crv.pan(vs_pan_dir)
    crv_list.append(crv)

vs_surf = RuledSurf(crv_list[0], crv_list[1])
gec.add_entity(vs_surf.to_iges())

'''Horizontal Stablizer'''
hs_x_start = tail_x_start + hs_x_offset
hs_y_start = hs_y_offset
hs_pan_dir = (hs_x_start, hs_y_start, 0)

foil = ['M6', 'M6']
z_offset = np.array([0, 4.4])
length = np.array([2.1, 1.4])
sweep_back = np.array([0, 12.0])
twist = np.zeros(2)
dihedral = np.array([0, -4])
twist_pos = np.array([0.25, 0.25])
y_ref = np.zeros(2)
thickness_factor = np.ones(2) * 1.5

crv_list = []
for k in range(2):
    wp = WingProfile.from_geom_param(foil[k], z_offset[k], length[k], sweep_back[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness_factor[k])
    crv = wp.nurbs_rep()
    crv.pan(hs_pan_dir)
    crv_list.append(crv)

hs_surf = RuledSurf(crv_list[0], crv_list[1])
hsm = deepcopy(hs_surf)
hsm.mirror('Z')
gec.add_entity(hs_surf.to_iges())
gec.add_entity(hsm.to_iges())

gec.write()
if auto_view:
    view(fn)
