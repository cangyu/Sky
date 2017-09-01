import numpy as np
import multiprocessing
import os
import math
from copy import deepcopy
from src.nurbs.curve import ClampedNURBSCrv, Line, Arc
from src.nurbs.surface import ClampedNURBSSurf, RuledSurf
from src.aircraft.wing import Wing, WingProfile
from src.aircraft.frame import BWBFrame
from src.msh.spacing import chebshev_dist_multi, hyperbolic_tangent, uniform, single_exponential, double_exponential
from src.msh.fluent import XF_MSH, BCType
from src.msh.tfi import LinearTFI3D, LinearTFI2D
from src.opt.latin import LatinHyperCube

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True


def calc_blk(_idx, _grid, _dist):
    print("PID: {}".format(os.getpid()))
    _u, _v, _w = _dist
    _grid.calc_grid(_u, _v, _w)
    print('Block {} calculation done!'.format(_idx))


'''Wing Planform'''
c_root = 11
c_mid = 5.6
c_tip = c_mid / 2.6
b_mid = 6.9
b_tip = 21
alpha_mid = 32
alpha_tip = 28

'''Profile distribution'''
inner_sec_num = 6
outer_sec_num = 9
u_mid = b_mid / b_tip
u_dist = chebshev_dist_multi([0, u_mid, 1], [inner_sec_num, outer_sec_num])
sec_num = len(u_dist)

frm = BWBFrame(c_root, c_mid, c_tip, b_mid, b_tip, alpha_mid, alpha_tip)
print(frm)
frm.show(u_dist)

'''Profile details'''
foil = ['SC20414', 'SC20414', 'SC20612', 'SC20712', 'SC20710', 'SC20710', 'SC20710', 'SC21010', 'SC21010', 'SC21006', 'SC20706', 'SC20706', 'SC20606', 'SC20406']
z_offset = list(map(frm.z, u_dist))
length = list(map(lambda u: frm.chord_len(u), u_dist))
sweep_back = list(map(lambda u: math.degrees(math.atan2(frm.x_front(u), frm.z(u))), u_dist))
twist = np.zeros(sec_num)
dihedral = np.zeros(sec_num)
twist_pos = np.full(sec_num, 0.25)
y_ref = np.zeros(sec_num)
thickness_factor = np.ones(sec_num)

'''LHS'''
n = sec_num * 10
twist_range = np.array([np.linspace(-0.50, 6.50, n),
                        np.linspace(-0.50, 6.50, n),
                        np.linspace(-2.50, 2.20, n),
                        np.linspace(-1.20, 4.00, n),
                        np.linspace(-0.50, 5.00, n),
                        np.linspace(-0.50, 5.00, n),
                        np.linspace(-0.50, 5.00, n),
                        np.linspace(- 2.50, 3.00, n),
                        np.linspace(- 2.50, 3.00, n),
                        np.linspace(- 3.00, 1.00, n),
                        np.linspace(- 1.20, 3.00, n),
                        np.linspace(- 1.20, 3.00, n),
                        np.linspace(- 0.80, 2.50, n),
                        np.linspace(- 1.20, 4.50, n)])

lhc = LatinHyperCube(twist_range)
sp = lhc.sample()

'''Record sampling results'''
fsp = open('sample_result.txt', 'w')
for _spk, param in enumerate(sp):
    fsp.write("Case {}:\n".format(_spk))
    fsp.write("{:^10} {:^10}\n".format('Profile', 'Twist'))
    for l in range(sec_num):
        fsp.write("{:^10}{:^10.2f}\n".format(l, param[l]))
    fsp.write('\n')
fsp.close()

'''Generate grids'''
for _spk, param in enumerate(sp):
    print("Building case {} grid ...".format(_spk))
    twist = np.copy(param)
    wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)
    wg.gen_grid('HWB_Wing_{}.msh'.format(_spk))
