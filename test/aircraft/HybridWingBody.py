import numpy as np
import math
from src.aircraft.wing import Wing
from src.aircraft.frame import BWBFrame
from src.msh.spacing import chebshev_dist_multi

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

s = 220
mac = 6.8
c_mid = 6
b_mid = 4.5
b_tip = 14
alpha_mid = 45
alpha_tip = 28
frm = BWBFrame.from_area_mac(s, mac, c_mid, b_mid, b_tip, alpha_mid, alpha_tip)

print(frm)
frm.show()

inner_sec_num = 6
outer_sec_num = 8
u_mid = b_mid / b_tip
u_dist = chebshev_dist_multi([0, u_mid, 1], [inner_sec_num, outer_sec_num])
sec_num = len(u_dist)

foil = ['NACA0012', 'NACA0012', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6']
z_offset = list(map(frm.z, u_dist))
length = list(map(lambda u: frm.x_tail(u) - frm.x_front(u), u_dist))
sweep_back = list(map(lambda u: math.degrees(math.atan2(frm.x_front(u), frm.z(u))), u_dist))
twist = np.array([3, 2.5, 2.2, 2, 2, 1.8, 1.6, 1.1, 1, 1, 0.5, 0.2, 0])
dihedral = np.array([0, 1, 1.1, 1.2, 2, 2, 2, 2, 2, 2, 2.5, 3, 3])
twist_pos = np.full(sec_num, 0.25)
y_ref = np.zeros(sec_num)
thickness_factor = np.ones(sec_num)

wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)

fn = 'HWB.igs'
model = wg.iges_model(fn)
model.write()

if auto_view:
    view(fn)
