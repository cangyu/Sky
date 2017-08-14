import numpy as np
import math
from copy import deepcopy
from src.aircraft.wing import Wing, Airfoil, WingProfile
from src.aircraft.frame import BWBFrame, chebshev_dist_multi
from src.nurbs.utility import pnt_dist
from src.nurbs.curve import ClampedNURBSCrv, Line, Arc, Spline
from src.nurbs.surface import Coons, BilinearSurf, ClampedNURBSSurf, RuledSurf
from src.iges.iges_core import IGES_Model
from src.msh.tfi import LinearTFI3D
from src.msh.plot3d import PLOT3D_Block, PLOT3D
from src.msh.fluent import XF_MSH, BCType
from src.msh.spacing import single_exponential, double_exponential, hyperbolic_tangent, hyperbolic_sine

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
fsf = RuledSurf(crv_tip, crv_far)

brk_root = np.array([0.44, 0.56])
brk_tip = brk_root
brk_far = brk_root

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
     Line(p[29], p[35]),  # 49
     Line(p[26], p[28]),  # 50
     Line(p[27], p[29]),  # 51
     Line(p[30], p[32]),  # 52
     Line(p[31], p[33]),  # 53
     Line(p[32], p[34]),  # 54
     Line(p[33], p[35])  # 55
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

model.add_entity(wsf.to_iges())
model.add_entity(fsf.to_iges())
model.write()
if auto_view:
    view(fn)

n0 = 50
n1 = 50
n2 = 50
n3 = 50
n4 = 50
n5 = 50
n6 = 50
n7 = 50
n = np.array([n0, n1, n2, n3, n4, n5, n6, n7])

u0 = np.linspace(0, 1, n0)
u1 = np.linspace(0, 1, n1)
u2 = np.linspace(0, 1, n2)
u3 = np.linspace(0, 1, n3)
u4 = np.linspace(0, 1, n4)
u5 = np.linspace(0, 1, n5)
u6 = np.linspace(0, 1, n6)
u7 = np.linspace(0, 1, n7)
knot_dist = [u0, u1, u2, u3, u4, u5, u6, u7]

'''Construct blocks'''
b0_s1 = Coons(l[35], l[37], c[0], l[20])
b0_s2 = Coons(l[26], l[29], l[13], l[12])
b0_s3 = Coons(c[0], l[13], l[1], l[5])
b0_s4 = Coons(l[20], l[12], l[0], l[4])
b0_s5 = Coons(l[1], l[0], l[35], l[26])
b0_s6 = Coons(l[5], l[4], l[37], l[29])

b0_tfi_grid = LinearTFI3D(lambda v, w: b0_s1(v, w),
                          lambda v, w: b0_s2(v, w),
                          lambda w, u: b0_s3(w, u),
                          lambda w, u: b0_s4(w, u),
                          lambda u, v: b0_s5(u, v),
                          lambda u, v: b0_s6(u, v))

b1_s1 = Coons(l[39], l[40], c[1], c[0])
b1_s2 = Coons(l[27], l[30], l[14], l[13])
b1_s3 = Coons(c[1], l[14], l[2], l[6])
b1_s4 = Coons(c[0], l[13], l[1], l[5])
b1_s5 = Coons(l[2], l[1], l[39], l[27])
b1_s6 = Coons(l[6], l[5], l[40], l[30])

b1_tfi_grid = LinearTFI3D(lambda v, w: b1_s1(v, w),
                          lambda v, w: b1_s2(v, w),
                          lambda w, u: b1_s3(w, u),
                          lambda w, u: b1_s4(w, u),
                          lambda u, v: b1_s5(u, v),
                          lambda u, v: b1_s6(u, v))

b2_s1 = Coons(l[2], l[6], c[1], l[14])
b2_s2 = Coons(l[3], l[7], l[21], l[15])
b2_s3 = Coons(c[1], l[21], l[36], l[38])
b2_s4 = Coons(l[14], l[15], l[28], l[31])
b2_s5 = Coons(l[36], l[28], l[2], l[3])
b2_s6 = Coons(l[38], l[31], l[6], l[7])

b2_tfi_grid = LinearTFI3D(lambda v, w: b2_s1(v, w),
                          lambda v, w: b2_s2(v, w),
                          lambda w, u: b2_s3(w, u),
                          lambda w, u: b2_s4(w, u),
                          lambda u, v: b2_s5(u, v),
                          lambda u, v: b2_s6(u, v))

b3_s1 = Coons(l[37], l[41], l[23], l[22])
b3_s2 = Coons(l[29], l[32], l[17], l[16])
b3_s3 = Coons(l[23], l[17], l[5], l[9])
b3_s4 = Coons(l[22], l[16], l[4], l[8])
b3_s5 = Coons(l[5], l[4], l[37], l[29])
b3_s6 = Coons(l[9], l[8], l[41], l[32])

b3_tfi_grid = LinearTFI3D(lambda v, w: b3_s1(v, w),
                          lambda v, w: b3_s2(v, w),
                          lambda w, u: b3_s3(w, u),
                          lambda w, u: b3_s4(w, u),
                          lambda u, v: b3_s5(u, v),
                          lambda u, v: b3_s6(u, v))

b4_s1 = Coons(l[40], l[42], l[24], l[23])
b4_s2 = Coons(l[30], l[33], l[18], l[17])
b4_s3 = Coons(l[24], l[18], l[6], l[10])
b4_s4 = Coons(l[23], l[17], l[5], l[9])
b4_s5 = Coons(l[6], l[5], l[40], l[30])
b4_s6 = Coons(l[10], l[9], l[42], l[33])

b4_tfi_grid = LinearTFI3D(lambda v, w: b4_s1(v, w),
                          lambda v, w: b4_s2(v, w),
                          lambda w, u: b4_s3(w, u),
                          lambda w, u: b4_s4(w, u),
                          lambda u, v: b4_s5(u, v),
                          lambda u, v: b4_s6(u, v))

b5_s1 = Coons(l[6], l[10], l[24], l[18])
b5_s2 = Coons(l[7], l[11], l[25], l[19])
b5_s3 = Coons(l[24], l[25], l[38], l[43])
b5_s4 = Coons(l[18], l[19], l[31], l[34])
b5_s5 = Coons(l[38], l[31], l[6], l[7])
b5_s6 = Coons(l[43], l[34], l[10], l[11])

b5_tfi_grid = LinearTFI3D(lambda v, w: b5_s1(v, w),
                          lambda v, w: b5_s2(v, w),
                          lambda w, u: b5_s3(w, u),
                          lambda w, u: b5_s4(w, u),
                          lambda u, v: b5_s5(u, v),
                          lambda u, v: b5_s6(u, v))

ts1 = ClampedNURBSSurf.split(wsf, brk_root, [])
s0 = ts1[0][0]
s1 = ts1[1][0]
s2 = ts1[2][0]

ts2 = ClampedNURBSSurf.split(fsf, brk_tip, [])
s3 = ts2[0][0]
s4 = ts2[1][0]
s5 = ts2[2][0]

s = [s0, s1, s2, s3, s4, s5]

b6_s1 = Coons(l[35], l[37], c[0], l[20])
b6_s2 = Coons(l[26], l[29], l[13], l[12])
b6_s3 = Coons(c[0], l[13], l[1], l[5])
b6_s4 = Coons(l[20], l[12], l[0], l[4])
b6_s5 = Coons(l[1], l[0], l[35], l[26])
b6_s6 = Coons(l[5], l[4], l[37], l[29])

b6_tfi_grid = LinearTFI3D(lambda v, w: b6_s1(v, w),
                          lambda v, w: b6_s2(v, w),
                          lambda w, u: b6_s3(w, u),
                          lambda w, u: b6_s4(w, u),
                          lambda u, v: b6_s5(u, v),
                          lambda u, v: b6_s6(u, v))

b7_s1 = Coons(l[39], l[40], c[1], c[0])
b7_s2 = Coons(l[27], l[30], l[14], l[13])
b7_s3 = Coons(c[1], l[14], l[2], l[6])
b7_s4 = Coons(c[0], l[13], l[1], l[5])
b7_s5 = Coons(l[2], l[1], l[39], l[27])
b7_s6 = Coons(l[6], l[5], l[40], l[30])

b7_tfi_grid = LinearTFI3D(lambda v, w: b7_s1(v, w),
                          lambda v, w: b7_s2(v, w),
                          lambda w, u: b7_s3(w, u),
                          lambda w, u: b7_s4(w, u),
                          lambda u, v: b7_s5(u, v),
                          lambda u, v: b7_s6(u, v))

b8_s1 = Coons(l[2], l[6], c[1], l[14])
b8_s2 = Coons(l[3], l[7], l[21], l[15])
b8_s3 = Coons(c[1], l[21], l[36], l[38])
b8_s4 = Coons(l[14], l[15], l[28], l[31])
b8_s5 = Coons(l[36], l[28], l[2], l[3])
b8_s6 = Coons(l[38], l[31], l[6], l[7])

b8_tfi_grid = LinearTFI3D(lambda v, w: b8_s1(v, w),
                          lambda v, w: b8_s2(v, w),
                          lambda w, u: b8_s3(w, u),
                          lambda w, u: b8_s4(w, u),
                          lambda u, v: b8_s5(u, v),
                          lambda u, v: b8_s6(u, v))

b9_s1 = Coons(l[37], l[41], l[23], l[22])
b9_s2 = Coons(l[29], l[32], l[17], l[16])
b9_s3 = Coons(l[23], l[17], l[5], l[9])
b9_s4 = Coons(l[22], l[16], l[4], l[8])
b9_s5 = Coons(l[5], l[4], l[37], l[29])
b9_s6 = Coons(l[9], l[8], l[41], l[32])

b9_tfi_grid = LinearTFI3D(lambda v, w: b9_s1(v, w),
                          lambda v, w: b9_s2(v, w),
                          lambda w, u: b9_s3(w, u),
                          lambda w, u: b9_s4(w, u),
                          lambda u, v: b9_s5(u, v),
                          lambda u, v: b9_s6(u, v))

b10_s1 = Coons(l[40], l[42], l[24], l[23])
b10_s2 = Coons(l[30], l[33], l[18], l[17])
b10_s3 = Coons(l[24], l[18], l[6], l[10])
b10_s4 = Coons(l[23], l[17], l[5], l[9])
b10_s5 = Coons(l[6], l[5], l[40], l[30])
b10_s6 = Coons(l[10], l[9], l[42], l[33])

b10_tfi_grid = LinearTFI3D(lambda v, w: b10_s1(v, w),
                           lambda v, w: b10_s2(v, w),
                           lambda w, u: b10_s3(w, u),
                           lambda w, u: b10_s4(w, u),
                           lambda u, v: b10_s5(u, v),
                           lambda u, v: b10_s6(u, v))

b11_s1 = Coons(l[6], l[10], l[24], l[18])
b11_s2 = Coons(l[7], l[11], l[25], l[19])
b11_s3 = Coons(l[24], l[25], l[38], l[43])
b11_s4 = Coons(l[18], l[19], l[31], l[34])
b11_s5 = Coons(l[38], l[31], l[6], l[7])
b11_s6 = Coons(l[43], l[34], l[10], l[11])

b11_tfi_grid = LinearTFI3D(lambda v, w: b11_s1(v, w),
                           lambda v, w: b11_s2(v, w),
                           lambda w, u: b11_s3(w, u),
                           lambda w, u: b11_s4(w, u),
                           lambda u, v: b11_s5(u, v),
                           lambda u, v: b11_s6(u, v))
