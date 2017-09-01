import numpy as np
import multiprocessing
import os
import math
from src.aircraft.wing import Wing
from src.aircraft.frame import BWBFrame
from src.msh.spacing import chebshev_dist_multi
from src.opt.latin import LatinHyperCube

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


if __name__ == '__main__':
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
    n = sec_num // 2
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

    '''Generate wings'''
    wing_list = []
    for _idx, _param in enumerate(sp):
        twist = np.copy(_param)
        wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)
        wing_list.append(wg)

    '''Generate grids'''
    report("{} sets of grid to be generated.".format(n))
    core_num = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=core_num)
    for k, wing in enumerate(wing_list):
        pool.apply_async(build_case_grid, (k, wing))
    pool.close()
    pool.join()
    report("All grid generation work done!")
