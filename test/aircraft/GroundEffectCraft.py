import numpy as np
import math
from src.aircraft.wing import Wing
from src.nurbs.curve import Arc
from src.nurbs.surface import ExtrudedSurf

try:
    from src.misc.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True

len_inner = 4.5
len_outer = 4.2
cr = 3.3
sweep_back = math.radians(25)
dihedral = math.radians(1 + 2)


def x_front(u):
    return 0 if u <= 0.5 else len_outer * math.tan(sweep_back) * (u - 0.5) / 0.5


def x_tail(u):
    return cr


def y(u):
    return u * (len_inner + len_outer) * math.tan(dihedral)


def z(u):
    return u * (len_inner + len_outer)


foil = 'M6'
sec_num = 17
u_dist = np.linspace(0, 1, sec_num)
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

fn = "GroundEffectCraft_{}_{}.igs".format(sec_num, foil)
igs_model = wg.write(fn, mirror=True)

'''Body'''
body_radius = 1.45
body_len = 11.2
body_profile = Arc.from_2pnt((-3.2, body_radius, 0), (-3.2, -body_radius, 0), 180, (1, 0, 0))
body_half = ExtrudedSurf(body_profile, (body_len, 0, 0)).to_iges()
igs_model.add_entity(body_half)
body_profile_mirrow = Arc.from_2pnt((-3.2, body_radius, 0), (-3.2, -body_radius, 0), 180, (-1, 0, 0))
body_half_mirrow = ExtrudedSurf(body_profile_mirrow, (body_len, 0, 0)).to_iges()
igs_model.add_entity(body_half_mirrow)




igs_model.write()
if auto_view:
    view(fn)
