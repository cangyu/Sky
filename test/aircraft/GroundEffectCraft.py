import numpy as np
import math
from src.iges.iges_core import IGES_Model
from src.nurbs.curve import Arc, Line, ClampedNURBSCrv, ConicArc, Spline
from src.nurbs.surface import Coons, ExtrudedSurf
from src.aircraft.wing import Airfoil, WingProfile

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

fn = "GroundEffectCraft.igs"
gec = IGES_Model(fn)

fuselage_width = 2.93
fuselage_height = 3.09

hu = 2.472
hd = hu - fuselage_height
nl = 4.66
nfr = 0.05
r1 = 0.9
theta1 = 15
rtheta1 = math.radians(theta1)
l1 = 0.16

p0 = np.array([0, nfr, 0])
p1 = np.array([r1 * (1 - math.cos(rtheta1)), nfr + r1 * math.sin(rtheta1), 0])
t1 = np.array([math.cos(math.radians(90 - theta1)), math.sin(math.radians(90 - theta1)), 0])
p2 = np.array([nl - l1, hu, 0])
t2 = np.array([1.0, 0, 0])
p12 = np.array([1.8, 1.77, 0])
p3 = np.array([nl, hu, 0])

crv1 = Arc.from_2pnt(p0, p1, theta1, (0, 0, -1))
crv2 = ConicArc(p1, t1, p2, t2, p12)
crv3 = Line(p2, p3)
crv4 = ClampedNURBSCrv.joint([crv1, crv2, crv3], 3)
gec.add_entity(crv4.to_iges())

r3 = 0.3
theta2 = 62
rtheta2 = math.radians(theta2)
l2 = 0.36

p4 = np.array([0, -nfr, 0])
p5 = np.array([r3 * (1 - math.cos(rtheta2)), -(nfr + r3 * math.sin(rtheta2)), 0])
t5 = np.array([math.sin(rtheta2), -math.cos(rtheta2), 0])
p6 = np.array([nl - l2, hd, 0])
t6 = np.array([1.0, 0, 0])
p13 = np.array([2, -0.56, 0])
p7 = np.array([nl, hd, 0])

crv5 = Arc.from_2pnt(p4, p5, theta2, (0, 0, 1))
crv6 = ConicArc(p5, t5, p6, t6, p13)
crv7 = Line(p6, p7)
crv8 = ClampedNURBSCrv.joint([crv5, crv6, crv7], 3)
gec.add_entity(crv8.to_iges())

crv9 = Arc.from_2pnt(p0, p4, 180, (1, 0, 0))
gec.add_entity(crv9.to_iges())

p8 = np.array([0, 0, nfr])
p9 = np.array([nl, 0.5 * (p3[1] + p7[1]), fuselage_width / 2])

crv9l0 = Arc.from_2pnt(p0, p8, 90, (1, 0, 0))
crv9l1 = Arc.from_2pnt(p8, p4, 90, (1, 0, 0))


def nose_side_off(y):
    return nl, y, fuselage_width / 2 * math.sqrt(1 - math.pow((y - p9[1]) / (fuselage_height / 2), 2))


c10l0ydst = np.linspace(hu, p9[1], 40)
c10l1ydst = np.linspace(p9[1], hd, 40)
c10ydst = np.linspace(hu, hd, 100)

c10l0pnt = list(map(nose_side_off, c10l0ydst))
c10l1pnt = list(map(nose_side_off, c10l1ydst))
c10pnt = list(map(nose_side_off, c10ydst))

crv10l0 = Spline(np.copy(c10l0pnt))
crv10l1 = Spline(np.copy(c10l1pnt))
crv10 = Spline(np.copy(c10pnt))

crv11 = Arc.from_2pnt(p8, p9, 32, np.cross(p8 - p9, (0, 0, 1)))
gec.add_entity(crv11.to_iges())

nose_surf1 = Coons(crv4, crv11, crv9l0, crv10l0)
nose_surf2 = Coons(crv11, crv8, crv9l1, crv10l1)
gec.add_entity(nose_surf1.to_iges())
gec.add_entity(nose_surf2.to_iges())

'''Fuselage'''
fuselage_len = 11.95
fuselage_surf = ExtrudedSurf(crv10, (fuselage_len, 0, 0))
gec.add_entity(fuselage_surf.to_iges())

'''Tail'''


'''Wing'''
foil = ['M6', 'M6']
z_offset = np.array([0, , 9.2])
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


gec.write()
if auto_view:
    view(fn)
