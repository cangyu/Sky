import numpy as np
import multiprocessing
import os
import math
from src.aircraft.wing import Wing
from src.aircraft.frame import HWBFrame
from src.msh.spacing import chebshev_dist_multi
from src.opt.latin import LatinHyperCube
from src.iges.iges_core import IGES_Model

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
cr = 12
spn = 21
fl = np.array([1.1, 2.3, 2.9])
alpha = np.radians([60, 62, 40, 28])
tl = np.array([1.4, 2.6, 2.5])
beta = np.radians([-18, -37, -20, 20])
frm = HWBFrame(cr, spn, fl, alpha, tl, beta)

fc = frm.fc
tc = frm.tc

fn = 'HWB_Wing.igs'
model = IGES_Model(fn)
model.add_entity(fc.to_iges())
model.add_entity(tc.to_iges())
model.write()

if auto_view:
    view(fn)

'''Profile distribution'''
inner_sec_num = 4
fusion_sec_num = 3
outer_sec_num = 5
u_dist = chebshev_dist_multi([0, 0.2, 0.4, 1], [inner_sec_num, fusion_sec_num, outer_sec_num])
sec_num = len(u_dist)

'''Profile details'''
foil = ['SC20414', 'SC20414', 'SC20612',
        'SC20712', 'SC20710',
        'SC20710', 'SC20710', 'SC21010', 'SC21010', 'SC21006', 'SC20706', 'SC20706', 'SC20606', 'SC20406']
