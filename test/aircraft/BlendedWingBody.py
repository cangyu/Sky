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

c_root = 14
c_mid = 6.8
c_tip = 1.3
b_mid = 8.1
b_tip = 21
alpha_mid = 32
alpha_tip = 28

inner_sec_num = 6
outer_sec_num = 9
u_mid = b_mid / b_tip
u_dist = chebshev_dist_multi([0, u_mid, 1], [inner_sec_num, outer_sec_num])
sec_num = len(u_dist)

frm = BWBFrame(c_root, c_mid, c_tip, b_mid, b_tip, alpha_mid, alpha_tip)
print(frm)
frm.show(u_dist)

foil = ['SC20414', 'SC20414', 'SC20612', 'SC20712', 'SC20710', 'SC20710',
        'SC20710', 'SC21010', 'SC21010', 'SC21006', 'SC20706', 'SC20706', 'SC20606', 'SC20406']
z_offset = list(map(frm.z, u_dist))
length = list(map(lambda u: frm.chord_len(u), u_dist))
sweep_back = list(map(lambda u: math.degrees(math.atan2(frm.x_front(u), frm.z(u))), u_dist))
twist = np.zeros(sec_num)
dihedral = np.zeros(sec_num)
twist_pos = np.full(sec_num, 0.25)
y_ref = np.zeros(sec_num)
thickness_factor = np.ones(sec_num)

wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)

fn = 'BWB.igs'
model = wg.iges_model(fn)
model.write()

if auto_view:
    view(fn)
