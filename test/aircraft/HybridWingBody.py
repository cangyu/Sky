import math
import os

import numpy as np
from matplotlib import pylab as plt

from src.aircraft.frame import HWBFrame
from src.aircraft.wing import Wing
from src.iges import IGES_Model
from src.msh.spacing import chebshev_dist_multi

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True


def report(msg):
    print('Process {} : {}'.format(os.getpid(), msg))


def build_case_grid(k, wing):
    report("Building case {} grid ...".format(k))
    wing.gen_grid('HWB_Wing_{}.msh'.format(k))


'''Wing Planform'''
spn = 21
cr = 16.667
ct = 1.38
outer_taper_ratio = 2.5
fl = np.array([2.1, 4, 2.9, 0])
alpha = np.radians([32, 56, 37.5, 28])
tl = np.array([1.5, 3.2, 0, 0])
beta = np.radians([-26, -40, 0, 0])
frm = HWBFrame(spn, cr, ct, fl, alpha, tl, beta, outer_taper=outer_taper_ratio)
print(frm)

'''Profile distribution'''
u = [0]
tmp = 0
for l in fl:
    tmp += l
    u.append(tmp)
tmp = 0
for l in tl:
    tmp += l
    u.append(tmp)
u = np.unique(u) / spn
u_dist = chebshev_dist_multi([0, u[1], u[2], u[3], u[4], u[5], u[6]], [3, 2, 3, 3, 3, 5])
frm.show(u_dist)

'''Profile details'''
sec_num = len(u_dist)
front_swp = np.zeros(sec_num)
for i in range(sec_num):
    tg = frm.front_crv(u_dist[i], 1)
    front_swp[i] = math.degrees(math.atan2(tg[0], tg[2]))

tc_3d = np.array([22, 22, 21, 19, 16, 14, 12, 11, 10, 8, 8, 8, 8, 8], float)
cl_3d = np.array([0.08, 0.10, 0.14, 0.18, 0.27, 0.38, 0.42, 0.45, 0.48, 0.47, 0.44, 0.29, 0.1, 0])
cl_2d = 1.1 * np.copy(cl_3d) / math.cos(math.radians(28)) ** 2
print(cl_2d)

foil = ['NACA14122',
        'NACA14022',
        'NACA63(4)-221',
        'NACA63(3)-218',
        'NLF(1)-0416',
        'NLF(1)-0414F',
        'SC(2)-0612',
        'SC(2)-0610',
        'SC(2)-0710',
        'SC(2)-0710',
        'SC(2)-0610',
        'SC(2)-0410',
        'NACA64A210',
        'NACA64A010']
z_offset = u_dist * spn
length = list(map(lambda _u: frm.chord_len(_u), u_dist))
sweep_back = list(map(lambda _u: math.degrees(math.atan2(frm.x_front(_u), frm.z(_u))), u_dist))
twist = np.zeros(sec_num)
dihedral = np.zeros(sec_num)
twist_pos = np.ones(sec_num)
y_ref = np.zeros(sec_num)
thickness_factor = np.ones(sec_num)

'''Show distribution'''
f, ax_arr = plt.subplots(2, sharex=True)
ax_arr[0].plot(u_dist * spn, tc_3d / 100)
ax_arr[0].set_ylim([0, tc_3d[0] / 100 * 1.1])
ax_arr[0].set_title('t/c in span-wise')
ax_arr[1].plot(u_dist * spn, cl_3d)
ax_arr[1].set_title('Cl in span-wise')
plt.show()

'''Initial grid'''
wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)

fn = 'HWB_Wing.igs'
model = IGES_Model(fn)
model.add_entity(wg.surf().to_iges())
model.write()
# if auto_view:
#     view(fn)

# wg.gen_grid('HWB_Wing.msh')
