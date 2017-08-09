import numpy as np
import math
from copy import deepcopy
from src.aircraft.wing import Wing, Airfoil, WingProfile
from src.aircraft.frame import BWBFrame, chebshev_dist_multi
from src.nurbs.utility import pnt_dist
from src.nurbs.curve import ClampedNURBSCrv, Line, Arc, Spline
from src.nurbs.surface import Coons, BilinearSurf
from src.iges.iges_core import IGES_Model
from src.msh.tfi import LinearTFI3D
from src.msh.plot3d import PLOT3D_Block, PLOT3D
from src.msh.fluent import XF_MSH, BCType

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

'''Wire Frame'''
fn = 'wireframe.igs'
model = IGES_Model(fn)

'''Shape parameters'''
c_root = 10
c_mid = 6
c_tip = 1.1
b_mid = 4.5
b_tip = 13
alpha_mid = 35
alpha_tip = 30
frm = BWBFrame(c_root, c_mid, c_tip, b_mid, b_tip, alpha_mid, alpha_tip)

'''Create wing from shape parameters'''
inner_sec_num = 6
outer_sec_num = 8
u_mid = b_mid / b_tip
u_dist = chebshev_dist_multi([0, u_mid, 1], [inner_sec_num, outer_sec_num])
sec_num = len(u_dist)

foil = ['M6' for x in range(sec_num)]
foil[0] = foil[1] = 'NACA0012'
z_offset = np.copy(list(map(frm.z, u_dist)))
length = np.copy(list(map(lambda u: frm.x_tail(u) - frm.x_front(u), u_dist)))
sweep_back = np.copy(list(map(lambda u: math.degrees(math.atan2(frm.x_front(u), frm.z(u))), u_dist)))
twist = np.array([3, 2.5, 2.2, 2, 2, 1.8, 1.6, 1.1, 1, 1, 0.5, 0.2, 0])
dihedral = np.array([0, 1, 1.1, 1.2, 2, 2, 2, 2, 2, 2, 2.5, 3, 3])
twist_pos = np.full(sec_num, 0.25)
y_ref = np.zeros(sec_num)
thickness_factor = np.ones(sec_num)
wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)
wsf = wg.surf()

'''Grid parameters'''
la = length[0]
lt = 30 * la
r = 10 * la
ispn = z_offset[-1]
ospn = 20 * ispn

crv_root = wsf.extract('V', 0)
crv_tip = wsf.extract('V', 1)
crv_far = WingProfile.from_geom_param(foil[-1], ispn + ospn, length[0], 0, 0, 0, thickness_factor=3).nurbs_rep()

brk_root = [0.44, 0.56]
brk_tip = [0.44, 0.56]
brk_far = [0.44, 0.56]

obrk_root = [0.3, 0.72]
obrk_tip = [0.3, 0.72]
obrk_far = [0.3, 0.72]

'''Points, lines, curves'''
p = np.zeros((36, 3))

p[2] = wsf(0, 0)
p[4] = wsf(1, 0)

p[0] = p[2]
p[0][1] += r

p[6] = p[4]
p[6][1] -= r

p[1] = p[0]
p[1][0] += lt

p[3] = p[2]
p[3][0] += lt

p[5] = p[4]
p[5][0] += lt

p[7] = p[6]
p[7][0] += lt

p[9] = p[1]
p[9][2] += ispn

p[11] = p[3]
p[11][2] += ispn

p[13] = p[5]
p[13][2] += ispn

p[15] = p[7]
p[15][2] += ispn

p[8] = p[0]
p[8][2] += ispn

p[10] = wsf(0, 1)
p[12] = wsf(1, 1)

p[14] = p[6]
p[14][2] += ispn

p[16] = p[8]
p[16][2] += ospn

p[17] = p[9]
p[17][2] += ospn

p[22] = p[14]
p[22][2] += ospn

p[23] = p[15]
p[23][2] += ospn

p[18] = crv_far.start
p[20] = crv_far.end

p[19] = p[18]
p[19][0] += lt

p[21] = p[20]
p[21][0] += lt

p[24] = crv_root(brk_root[0])
p[25] = crv_root(brk_root[1])

p[26] = crv_tip(brk_tip[0])
p[27] = crv_tip(brk_tip[1])

p[28] = crv_far(brk_far[0])
p[29] = crv_far(brk_far[1])

outer_root = Arc.from_2pnt(p[0], p[6], 180, (0, 0, 1))
outer_tip = Arc.from_2pnt(p[8], p[14], 180, (0, 0, 1))
outer_far = Arc.from_2pnt(p[16], p[22], 180, (0, 0, 1))

p[30] = outer_root(obrk_root[0])
p[31] = outer_root(obrk_root[1])
p[32] = outer_tip(obrk_tip[0])
p[33] = outer_tip(obrk_tip[1])
p[34] = outer_far(obrk_far[0])
p[35] = outer_far(obrk_far[1])

l = [Line(p[0], p[1]),  # 0
     Line(p[2], p[3]),  # 1
     Line(p[4], p[5]),  # 2
     Line(p[6], p[7]),  # 3
     Line(p[8], p[9]),  # 4
     Line(p[10], p[11]),  # 5
     Line(p[12], p[13]),  # 6
     Line(p[14], p[15]),  # 7
     Line(p[16], p[17]),  # 8
     Line(p[18], p[19]),  # 9
     Line(p[20], p[21]),  # 10
     Line(p[22], p[23]),  # 11
     Line(p[1], p[9]),  # 12
     Line(p[3], p[11]),  # 13
     Line(p[5], p[13]),  # 14
     Line(p[7], p[15]),  # 15
     Line(p[9], p[17]),  # 16
     Line(p[11], p[19]),  # 17
     Line(p[13], p[21]),  # 18
     Line(p[15], p[23]),  # 19
     Line(p[0], p[8]),  # 20
     Line(p[6], p[14]),  # 21
     Line(p[8], p[16]),  # 22
     Line(p[10], p[18]),  # 23
     Line(p[12], p[20]),  # 24
     Line(p[14], p[22]),  # 25
     Line(p[3], p[1]),  # 26
     Line(p[5], p[3]),  # 27
     Line(p[5], p[7]),  # 28
     Line(p[11], p[9]),  # 29
     Line(p[13], p[11]),  # 30
     Line(p[13], p[15]),  # 31
     Line(p[19], p[17]),  # 32
     Line(p[21], p[19]),  # 33
     Line(p[21], p[23]),  # 34
     Line(p[2], p[0]),  # 35
     Line(p[4], p[6]),  # 36
     Line(p[10], p[8]),  # 37
     Line(p[12], p[14]),  # 38
     Line(p[4], p[2]),  # 39
     Line(p[12], p[10]),  # 40
     Line(p[18], p[16]),  # 41
     Line(p[20], p[18]),  # 42
     Line(p[20], p[22]),  # 43
     Line(p[24], p[30]),  # 44
     Line(p[25], p[31]),  # 45
     Line(p[26], p[32]),  # 46
     Line(p[27], p[33]),  # 47
     Line(p[28], p[34]),  # 48
     Line(p[29], p[35])  # 49
     ]

c = [wsf.extract('U', 0), wsf.extract('U', 1)]
c2, c3, c4 = ClampedNURBSCrv.split(outer_root, obrk_root)
c5, c6, c7 = ClampedNURBSCrv.split(outer_tip, obrk_tip)
c8, c9, c10 = ClampedNURBSCrv.split(outer_far, obrk_far)
c11, c12, c13 = ClampedNURBSCrv.split(crv_root, brk_root)
c14, c15, c16 = ClampedNURBSCrv.split(crv_tip, brk_tip)
c17, c18, c19 = ClampedNURBSCrv.split(crv_far, brk_far)

c4.reverse()
c7.reverse()
c10.reverse()
c13.reverse()
c16.reverse()
c19.reverse()

c.append(c2)
c.append(c3)
c.append(c4)

c.append(c5)
c.append(c6)
c.append(c7)

c.append(c8)
c.append(c9)
c.append(c10)

c.append(c11)
c.append(c12)
c.append(c13)

c.append(c14)
c.append(c15)
c.append(c16)

c.append(c17)
c.append(c18)
c.append(c19)

for _ln in l:
    model.add_entity(_ln.to_iges())

for _cv in c:
    model.add_entity(_cv.to_iges())

model.write()
if auto_view:
    view(fn)
