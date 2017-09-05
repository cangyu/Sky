import numpy as np
import multiprocessing
import os
import math
from src.aircraft.wing import Wing
from src.aircraft.frame import HWBFrame
from src.msh.spacing import chebshev_dist_multi
from src.opt.latin import LatinHyperCube
from src.iges.iges_core import IGES_Model
from matplotlib import pylab as plt

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
print("Outer Span: {} m".format(fl[-1]))
print("Outer AR: {}".format(2 * fl[-1] / (ct * (1 + outer_taper_ratio) / 2)))

# fn = 'HWB_Wing.igs'
# model = IGES_Model(fn)
# model.add_entity(frm.front_crv.to_iges())
# model.add_entity(frm.tail_crv.to_iges())
# model.add_entity(frm.root_line.to_iges())
# model.add_entity(frm.tip_line.to_iges())
# model.add_entity(frm.planform_surf.to_iges())
# model.write()
# if auto_view:
#     view(fn)

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
swp = np.array([0, 15, 21, 32, 56, 56, 56], float)
tc = np.array([24, 24, 23, 22, 18, 14, 13, 12, 11, 10, 9, 8, 8, 8], float)
cl = np.array([0.08, 0.10, 0.14, 0.18, 0.27, 0.38, 0.42, 0.45, 0.48, 0.47, 0.42, 0.28, 0.09, 0], float)
print(1.1 * cl / math.cos(math.radians(28)) ** 2)

'''Show distribution'''
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(u_dist * spn, tc / 100)
axarr[0].set_ylim([0, tc[0] / 100 * 1.1])
axarr[0].set_title('t/c in span-wise')
axarr[1].plot(u_dist * spn, cl)
axarr[1].set_title('Cl in span-wise')
plt.show()

foil = ['NACA23024', 'SC20414',
        'SC20612',
        'SC20712', 'SC20710',
        'SC20710', 'SC20710',
        'SC21006', 'SC20706', 'SC20706',
        'SC20606', 'SC20406', 'SC21010', 'SC21006', 'SC21006']
