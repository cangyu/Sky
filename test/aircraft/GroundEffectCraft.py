import numpy as np
import math
from src.aircraft.wing import Wing

try:
    from src.misc.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True

len_inner = 2
len_outer = 2
cr = 3.3
sweep_back = math.radians(20)
dihedral = math.radians(1 + 2)


def x_front(u):
    return 0 if u <= 0.5 else len_outer * math.tan(sweep_back) * (u - 0.5) / 0.5


def x_tail():
    return cr


def y(u):
    return u * (len_inner + len_outer) * math.tan(dihedral)


def z(u):
    return u * (len_inner + len_outer)

foil = 'M6'
u_dist = np.array([0, 0.1, 0.3, 0.5, 0.6, 0.8, 1.0])
sec_num = len(u_dist)
airfoil_list = []
thickness = np.ones(sec_num, float)
z_dist = np.empty(sec_num, float)
xf_dist = np.empty(sec_num, float)
yf_dist = np.empty(sec_num, float)
xt_dist = np.empty(sec_num, float)
yt_dist = np.empty(sec_num, float)

for i in range(sec_num):
    airfoil_list.append(foil)
    z_dist[i] = z(u_dist[i])
    xf_dist[i] = x_front(u_dist[i])
    xt_dist[i] = x_tail(u_dist[i])
    yf_dist[i] = yt_dist[i] = y(u_dist[i])

wg = Wing(airfoil_list, thickness, z_dist, xf_dist, yf_dist, xt_dist, yt_dist)

fn = "GroundEffectCraft_{}_{}.igs".format(n, foil)
wg.write(fn, mirror=True)

if auto_view:
    view(fn)
