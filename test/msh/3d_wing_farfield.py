import math
from copy import deepcopy
from src.aircraft.wing import Wing, WingProfile
from src.aircraft.frame import BWBFrame
from src.nurbs.curve import ClampedNURBSCrv, Line, Arc
from src.nurbs.surface import ClampedNURBSSurf, RuledSurf
from src.iges.iges_core import IGES_Model
from src.msh.tfi import LinearTFI3D, LinearTFI2D
from src.msh.plot3d import PLOT3D_Block, PLOT3D
from src.msh.fluent import XF_MSH, BCType
from src.msh.spacing import *

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True


def report_process(idx):
    print('Block {} calculation done!'.format(idx))


'''Plot3D Representation'''
p3d_grid = PLOT3D()

'''Shape parameters'''
c_root = 14
c_mid = 6.8
c_tip = 1.3
b_mid = 8.1
b_tip = 21
alpha_mid = 32
alpha_tip = 28
frm = BWBFrame(c_root, c_mid, c_tip, b_mid, b_tip, alpha_mid, alpha_tip)
print(frm)

inner_sec_num = 6
outer_sec_num = 9
u_mid = b_mid / b_tip
u_dist = chebshev_dist_multi([0, u_mid, 1], [inner_sec_num, outer_sec_num])

sec_num = len(u_dist)
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
wsf = wg.surf()

'''Grid parameters'''
la = length[0]
lt = 30 * la
r = 10 * la
inner_spn = z_offset[-1]
outer_spn = 20 * inner_spn

crv_root = wsf.extract('V', 0)
crv_tip = wsf.extract('V', 1)
crv_far = WingProfile.from_geom_param(foil[-1], inner_spn + outer_spn, length[0], 0, 0, 0, thickness_factor=3).nurbs_rep()
fsf = RuledSurf(crv_tip, crv_far)

brk_root = np.array([0.44, 0.56])
brk_tip = brk_root
brk_far = brk_root

obrk_root = [0.3, 0.72]
obrk_tip = [0.3, 0.72]
obrk_far = [0.3, 0.72]

'''Points, lines, curves, surfs'''
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
p[9][2] += inner_spn

p[11] = p[3]
p[11][2] += inner_spn

p[13] = p[5]
p[13][2] += inner_spn

p[15] = p[7]
p[15][2] += inner_spn

p[8] = p[0]
p[8][2] += inner_spn

p[10] = wsf(0, 1)
p[12] = wsf(1, 1)

p[14] = p[6]
p[14][2] += inner_spn

p[16] = p[8]
p[16][2] += outer_spn

p[17] = p[9]
p[17][2] += outer_spn

p[22] = p[14]
p[22][2] += outer_spn

p[23] = p[15]
p[23][2] += outer_spn

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
c20 = wsf.extract('U', brk_root[0])
c21 = wsf.extract('U', brk_root[1])
c22 = fsf.extract('U', brk_tip[0])
c23 = fsf.extract('U', brk_tip[1])

c3.reverse()
c4.reverse()
c6.reverse()
c7.reverse()
c9.reverse()
c10.reverse()
c12.reverse()
c13.reverse()
c15.reverse()
c16.reverse()
c18.reverse()
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

c.append(c20)
c.append(c21)

c.append(c22)
c.append(c23)

ts1 = ClampedNURBSSurf.split(wsf, brk_root, [])
s0 = ts1[0][0]
s1 = ts1[1][0]
s2 = ts1[2][0]

ts2 = ClampedNURBSSurf.split(fsf, brk_tip, [])
s3 = ts2[0][0]
s4 = ts2[1][0]
s5 = ts2[2][0]

s = [s0, s1, s2, s3, s4, s5]

'''Wire Frame'''
fn = 'WireFrame.igs'
model = IGES_Model(fn)

for _ln in l:
    model.add_entity(_ln.to_iges())

for _cv in c:
    model.add_entity(_cv.to_iges())

for _sf in s:
    model.add_entity(_sf.to_iges())

model.write()
# if auto_view:
#     view(fn)

n0 = 60
n1 = 60
n2 = 60
n3 = 60
n4 = 60
n5 = 60
n6 = 60
n7 = 60

n = np.array([n0, n1, n2, n3, n4, n5, n6, n7])

u0 = hyperbolic_tangent(n0, 8)
u1 = double_exponential(n1, 0.5, 1.5, 0.5)
u2 = uniform(n2)
u3 = single_exponential(n3, 5)
u4 = hyperbolic_tangent(n4, 5)
u5 = double_exponential(n5, 0.5, 1.2, 0.5)
u6 = double_exponential(n6, 0.5, 1.5, 0.5)
u7 = uniform(n7)

knot_dist = [u0, u1, u2, u3, u4, u5, u6, u7]

'''Construct blocks'''
b0_tfi_grid = LinearTFI3D.from_edges(l[1], l[26], l[0], l[35], l[5], l[29], l[4], l[37], c[0], l[13], l[12], l[20])
b0_tfi_grid.calc_grid(knot_dist[3], knot_dist[0], knot_dist[2])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b0_tfi_grid.get_grid()))
report_process(0)

b1_tfi_grid = LinearTFI3D.from_edges(l[2], l[27], l[1], l[39], l[6], l[30], l[5], l[40], c[1], l[14], l[13], c[0])
b1_tfi_grid.calc_grid(knot_dist[3], knot_dist[7], knot_dist[2])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b1_tfi_grid.get_grid()))
report_process(1)

b2_tfi_grid = LinearTFI3D.from_edges(l[36], l[3], l[28], l[2], l[38], l[7], l[31], l[6], c[1], l[21], l[15], l[14])
b2_tfi_grid.calc_grid(knot_dist[0], knot_dist[3], knot_dist[2])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b2_tfi_grid.get_grid()))
report_process(2)

b3_tfi_grid = LinearTFI3D.from_edges(l[5], l[29], l[4], l[37], l[9], l[32], l[8], l[41], l[23], l[17], l[16], l[22])
b3_tfi_grid.calc_grid(knot_dist[3], knot_dist[0], knot_dist[4])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b3_tfi_grid.get_grid()))
report_process(3)

b4_tfi_grid = LinearTFI3D.from_edges(l[6], l[30], l[5], l[40], l[10], l[33], l[9], l[42], l[24], l[18], l[17], l[23])
b4_tfi_grid.calc_grid(knot_dist[3], knot_dist[7], knot_dist[4])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b4_tfi_grid.get_grid()))
report_process(4)

b5_tfi_grid = LinearTFI3D.from_edges(l[38], l[7], l[31], l[6], l[43], l[11], l[34], l[10], l[24], l[25], l[19], l[18])
b5_tfi_grid.calc_grid(knot_dist[0], knot_dist[3], knot_dist[4])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b5_tfi_grid.get_grid()))
report_process(5)

b6_s1 = deepcopy(s[0])
b6_s2 = LinearTFI2D(c[2], l[20], c[5], l[52])
b6_s3 = LinearTFI2D(c[0], l[35], l[20], l[37])
b6_s4 = LinearTFI2D(c[20], l[44], l[52], l[46])
b6_s5 = LinearTFI2D(l[35], c[11], l[44], c[2])
b6_s6 = LinearTFI2D(l[37], c[14], l[46], c[5])

b6_tfi_grid = LinearTFI3D(lambda v, w: b6_s1(v, w),
                          lambda v, w: b6_s2(v, w),
                          lambda w, u: b6_s3(w, u),
                          lambda w, u: b6_s4(w, u),
                          lambda u, v: b6_s5(u, v),
                          lambda u, v: b6_s6(u, v))

b6_tfi_grid.calc_grid(knot_dist[0], knot_dist[1], knot_dist[2])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b6_tfi_grid.get_grid()))
report_process(6)

b7_s1 = LinearTFI2D(l[45], c[21], l[47], l[53])
b7_s2 = LinearTFI2D(l[44], c[20], l[46], l[52])
b7_s3 = deepcopy(s[1])
b7_s3.reverse('U')
b7_s3.swap()
b7_s4 = LinearTFI2D(l[53], c[3], l[52], c[6])
b7_s5 = LinearTFI2D(c[12], l[45], c[3], l[44])
b7_s6 = LinearTFI2D(c[15], l[47], c[6], l[46])

b7_tfi_grid = LinearTFI3D(lambda v, w: b7_s1(v, w),
                          lambda v, w: b7_s2(v, w),
                          lambda w, u: b7_s3(w, u),
                          lambda w, u: b7_s4(w, u),
                          lambda u, v: b7_s5(u, v),
                          lambda u, v: b7_s6(u, v))

b7_tfi_grid.calc_grid(knot_dist[5], knot_dist[0], knot_dist[2])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b7_tfi_grid.get_grid()))
report_process(7)

b8_s1 = LinearTFI2D(l[36], c[1], l[38], l[21])
b8_s2 = LinearTFI2D(l[45], c[21], l[47], l[53])
b8_s3 = deepcopy(s[2])
b8_s3.reverse('U')
b8_s3.swap()
b8_s4 = LinearTFI2D(l[21], c[4], l[53], c[7])
b8_s5 = LinearTFI2D(c[13], l[36], c[4], l[45])
b8_s6 = LinearTFI2D(c[16], l[38], c[7], l[47])

b8_tfi_grid = LinearTFI3D(lambda v, w: b8_s1(v, w),
                          lambda v, w: b8_s2(v, w),
                          lambda w, u: b8_s3(w, u),
                          lambda w, u: b8_s4(w, u),
                          lambda u, v: b8_s5(u, v),
                          lambda u, v: b8_s6(u, v))

b8_tfi_grid.calc_grid(knot_dist[6], knot_dist[0], knot_dist[2])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b8_tfi_grid.get_grid()))
report_process(8)

b9_s1 = deepcopy(s[3])
b9_s2 = LinearTFI2D(c[5], l[22], c[8], l[54])
b9_s3 = LinearTFI2D(l[23], l[37], l[22], l[41])
b9_s4 = LinearTFI2D(c[22], l[46], l[54], l[48])
b9_s5 = LinearTFI2D(l[37], c[14], l[46], c[5])
b9_s6 = LinearTFI2D(l[41], c[17], l[48], c[8])

b9_tfi_grid = LinearTFI3D(lambda v, w: b9_s1(v, w),
                          lambda v, w: b9_s2(v, w),
                          lambda w, u: b9_s3(w, u),
                          lambda w, u: b9_s4(w, u),
                          lambda u, v: b9_s5(u, v),
                          lambda u, v: b9_s6(u, v))

b9_tfi_grid.calc_grid(knot_dist[0], knot_dist[1], knot_dist[4])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b9_tfi_grid.get_grid()))
report_process(9)

b10_s1 = LinearTFI2D(l[47], c[23], l[49], l[55])
b10_s2 = LinearTFI2D(l[46], c[22], l[48], l[54])
b10_s3 = deepcopy(s[4])
b10_s3.reverse('U')
b10_s3.swap()
b10_s4 = LinearTFI2D(l[55], c[6], l[54], c[9])
b10_s5 = LinearTFI2D(c[15], l[47], c[6], l[46])
b10_s6 = LinearTFI2D(c[18], l[49], c[9], l[48])

b10_tfi_grid = LinearTFI3D(lambda v, w: b10_s1(v, w),
                           lambda v, w: b10_s2(v, w),
                           lambda w, u: b10_s3(w, u),
                           lambda w, u: b10_s4(w, u),
                           lambda u, v: b10_s5(u, v),
                           lambda u, v: b10_s6(u, v))

b10_tfi_grid.calc_grid(knot_dist[5], knot_dist[0], knot_dist[4])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b10_tfi_grid.get_grid()))
report_process(10)

b11_s1 = LinearTFI2D(l[38], l[24], l[43], l[25])
b11_s2 = LinearTFI2D(l[47], c[23], l[49], l[55])
b11_s3 = deepcopy(s[5])
b11_s3.reverse('U')
b11_s3.swap()
b11_s4 = LinearTFI2D(l[25], c[7], l[55], c[10])
b11_s5 = LinearTFI2D(c[16], l[38], c[7], l[47])
b11_s6 = LinearTFI2D(c[19], l[43], c[10], l[49])

b11_tfi_grid = LinearTFI3D(lambda v, w: b11_s1(v, w),
                           lambda v, w: b11_s2(v, w),
                           lambda w, u: b11_s3(w, u),
                           lambda w, u: b11_s4(w, u),
                           lambda u, v: b11_s5(u, v),
                           lambda u, v: b11_s6(u, v))

b11_tfi_grid.calc_grid(knot_dist[6], knot_dist[0], knot_dist[4])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b11_tfi_grid.get_grid()))
report_process(11)

b12_s1 = deepcopy(s[5])
b12_s1.reverse('U')
b12_s2 = deepcopy(s[3])
b12_s3 = LinearTFI2D(l[24], l[40], l[23], l[42])
b12_s4 = deepcopy(s[4])
b12_s4.reverse('U')
b12_s4.swap()
b12_s5 = LinearTFI2D(l[40], c[16], c[15], c[14])
b12_s6 = LinearTFI2D(l[42], c[19], c[18], c[17])

b12_tfi_grid = LinearTFI3D(lambda v, w: b12_s1(v, w),
                           lambda v, w: b12_s2(v, w),
                           lambda w, u: b12_s3(w, u),
                           lambda w, u: b12_s4(w, u),
                           lambda u, v: b12_s5(u, v),
                           lambda u, v: b12_s6(u, v))

b12_tfi_grid.calc_grid(knot_dist[7], knot_dist[6], knot_dist[4])
p3d_grid.add_block(PLOT3D_Block.build_from_3d(b12_tfi_grid.get_grid()))
report_process(12)

'''网格, 边界条件, 邻接关系'''
blk = [b0_tfi_grid.get_grid(),
       b1_tfi_grid.get_grid(),
       b2_tfi_grid.get_grid(),
       b3_tfi_grid.get_grid(),
       b4_tfi_grid.get_grid(),
       b5_tfi_grid.get_grid(),
       b6_tfi_grid.get_grid(),
       b7_tfi_grid.get_grid(),
       b8_tfi_grid.get_grid(),
       b9_tfi_grid.get_grid(),
       b10_tfi_grid.get_grid(),
       b11_tfi_grid.get_grid(),
       b12_tfi_grid.get_grid()]

bc = [(BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b0
      (BCType.Wall, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Symmetry, BCType.Interior),  # b1
      (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b2
      (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b3
      (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField),  # b4
      (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b5
      (BCType.Wall, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Symmetry, BCType.Interior),  # b6
      (BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b7
      (BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b8
      (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField),  # b9
      (BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b10
      (BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b11
      (BCType.Interior, BCType.Interior, BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField)  # b12
      ]

adj = [((6, 3), (0, 1), 1, True),
       ((0, 2), (0, 0), 0, False),
       ((1, 4), (0, 3), 1, False),
       ((0, 4), (0, 0), 0, False),
       ((0, 0), (0, 5), 1, False),
       ((0, 6), (3, 5), 0, False),

       ((0, 0), (1, 1), 1, False),
       ((1, 2), (0, 0), 0, False),
       ((2, 1), (1, 3), 1, True),
       ((0, 0), (1, 5), 1, False),
       ((1, 6), (4, 5), 0, False),

       ((2, 2), (0, 0), 0, False),
       ((8, 1), (2, 3), 1, True),
       ((2, 4), (0, 0), 0, False),
       ((0, 0), (2, 5), 1, False),
       ((2, 6), (5, 5), 0, False),

       ((9, 3), (3, 1), 1, True),
       ((3, 2), (0, 0), 0, False),
       ((4, 4), (3, 3), 1, False),
       ((3, 4), (0, 0), 0, False),
       ((3, 6), (0, 0), 0, False),

       ((12, 3), (4, 1), 1, True),
       ((4, 2), (0, 0), 0, False),
       ((5, 1), (4, 3), 1, True),
       ((4, 6), (0, 0), 0, False),

       ((5, 2), (0, 0), 0, False),
       ((11, 1), (5, 3), 1, True),
       ((5, 4), (0, 0), 0, False),
       ((5, 6), (0, 0), 0, False),

       ((0, 0), (6, 1), 1, False),
       ((6, 2), (0, 0), 0, False),
       ((6, 4), (7, 2), 0, True),
       ((0, 0), (6, 5), 1, False),
       ((6, 6), (9, 5), 0, False),

       ((8, 2), (7, 1), 1, False),
       ((0, 0), (7, 3), 1, False),
       ((7, 4), (0, 0), 0, False),
       ((0, 0), (7, 5), 1, False),
       ((7, 6), (10, 5), 0, False),

       ((0, 0), (8, 3), 1, False),
       ((8, 4), (0, 0), 0, False),
       ((0, 0), (8, 5), 1, False),
       ((8, 6), (11, 5), 0, False),

       ((12, 2), (9, 1), 1, False),
       ((9, 2), (0, 0), 0, False),
       ((9, 4), (10, 2), 0, True),
       ((9, 6), (0, 0), 0, False),

       ((11, 2), (10, 1), 1, False),
       ((12, 4), (10, 3), 1, False),
       ((10, 4), (0, 0), 0, False),
       ((10, 6), (0, 0), 0, False),

       ((12, 1), (11, 3), 1, True),
       ((11, 4), (0, 0), 0, False),
       ((11, 6), (0, 0), 0, False),

       ((0, 0), (12, 5), 1, False),
       ((12, 6), (0, 0), 0, False)]

'''构建MSH文件'''
p3d_grid.write('HWB_Wing.xyz')
msh = XF_MSH.from_str3d_multi(blk, bc, adj)
msh.save('HWB_Wing.msh')
