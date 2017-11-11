import math
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from nurbs import to_homogeneous, to_cartesian, Crv, ConicArc, Arc, Line, Coons, RuledSurf
from iges import Model
from wing import WingProfile
from misc import sqrt2


'''Wing Calculation'''
'''Area and ratio'''
s = 38.88
sir = 0.68
sor = 1.0 - sir
si = sir * s
so = sor * s

'''inner'''
icr = 5.8
ispn = si / icr
ieta = ispn / icr

'''outer'''
ocr = 3.9
oeta = 1.15
ospn = oeta * ocr
octp = 2 * so / ospn - ocr
ofswp = 22

'''Geom param'''
w1 = icr
h1 = ispn
w2 = 1.2
w3 = ocr
h2 = ospn
theta = math.radians(ofswp)
w4 = octp

'''Points'''
p0 = (0, 0)
p1 = (0, w1)
p2 = (h1, w1)
p3 = (h1, 0)
p4 = (h1, w2 + w3)
p5 = (h1, w2)
p6 = (h1 + h2, p4[1] - h2 * math.tan(theta))
p7 = (p6[0], p6[1] - w4)

'''Report'''
print('单侧机翼总面积：{} m^2'.format(s))
print('两段面积比：{} / {}'.format(sir, sor))
print('内翼')
print('面积：{} m^2'.format(si))
print('根弦长：{} m'.format(icr))
print('展长：{} m'.format(ispn))
print('展弦比：{}'.format(ieta))
print('外翼')
print('面积：{} m^2'.format(so))
print('根弦长：{} m'.format(ocr))
print('梢弦长：{} m'.format(octp))
print('展长：{} m'.format(ospn))
print('展弦比：{}'.format(oeta))
print('根梢比：{}'.format(ocr / octp))
print('前缘后掠角：{}'.format(ofswp))

'''Plot'''
figure, ax = plt.subplots()
ax.set_xlim(left=0, right=11)
ax.set_ylim(bottom=-2, top=7)

'''Line endings'''
line1 = [p0, p1]
line2 = [p1, p2]
line3 = [p2, p3]
line4 = [p3, p0]
line5 = [p4, p6]
line6 = [p5, p7]
line7 = [p6, p7]

(line1_xs, line1_ys) = zip(*line1)
(line2_xs, line2_ys) = zip(*line2)
(line3_xs, line3_ys) = zip(*line3)
(line4_xs, line4_ys) = zip(*line4)
(line5_xs, line5_ys) = zip(*line5)
(line6_xs, line6_ys) = zip(*line6)
(line7_xs, line7_ys) = zip(*line7)

'''Create the lines'''
ax.add_line(Line2D(line1_xs, line1_ys, linewidth=1, color='blue'))
ax.add_line(Line2D(line2_xs, line2_ys, linewidth=1, color='red'))
ax.add_line(Line2D(line3_xs, line3_ys, linewidth=1, color='blue'))
ax.add_line(Line2D(line4_xs, line4_ys, linewidth=1, color='red'))
ax.add_line(Line2D(line5_xs, line5_ys, linewidth=1, color='blue'))
ax.add_line(Line2D(line6_xs, line6_ys, linewidth=1, color='red'))
ax.add_line(Line2D(line7_xs, line7_ys, linewidth=1, color='blue'))

'''Show'''
plt.plot()
plt.show()

'''Vertical-Tail calculation'''
'''Design variable'''
s = 6.866
eta = 1.3
nbla = 1.4
fswp = 40

spn = math.sqrt(s * eta)
cr = (2 * s / spn) / (1 + nbla) * nbla
ct = cr / nbla

'''Geom param'''
w1 = cr
w2 = ct
h1 = spn
theta = math.radians(fswp)

'''Report'''
print("垂尾面积：{} m^2".format(s))
print("展弦比：{}".format(eta))
print("根梢比：{}".format(nbla))
print("根弦长：{} m".format(cr))
print("梢弦长：{} m".format(ct))
print("高度：{} m".format(spn))
print("后掠角：{}".format(fswp))

'''Points'''
p0 = (0, 0)
p1 = (w1, 0)
p3 = (h1 * math.tan(theta), h1)
p2 = (p3[0] + w2, p3[1])

'''Plot'''
figure, ax = plt.subplots()
# ax.set_xlim(left=-11, right=11)
# ax.set_ylim(bottom=-2, top=7)

'''Line endings'''
line1 = [p0, p1]
line2 = [p1, p2]
line3 = [p2, p3]
line4 = [p3, p0]

(line1_xs, line1_ys) = zip(*line1)
(line2_xs, line2_ys) = zip(*line2)
(line3_xs, line3_ys) = zip(*line3)
(line4_xs, line4_ys) = zip(*line4)

'''Create the lines'''
ax.add_line(Line2D(line1_xs, line1_ys, linewidth=1, color='blue'))
ax.add_line(Line2D(line2_xs, line2_ys, linewidth=1, color='red'))
ax.add_line(Line2D(line3_xs, line3_ys, linewidth=1, color='blue'))
ax.add_line(Line2D(line4_xs, line4_ys, linewidth=1, color='red'))

'''Show'''
plt.plot()
plt.gca().set_aspect('equal')
plt.show()

'''Design variable'''
s = 20.11
eta = 3.8
nbla = 1.7
fswp = 10

spn = math.sqrt(s * eta)
cr = (2 * s / spn) / (1 + nbla) * nbla
ct = cr / nbla

'''Geom param'''
w1 = cr
w2 = ct
h1 = spn / 2
theta = math.radians(fswp)

'''Report'''
print("平尾面积：{} m^2".format(s))
print("展弦比：{}".format(eta))
print("根梢比：{}".format(nbla))
print("根弦长：{} m".format(cr))
print("梢弦长：{} m".format(ct))
print("展长：{} m".format(spn))
print("后掠角：{}".format(fswp))

'''Points'''
p0 = (0, 0)
p1 = (0, w1)
p2 = (h1, w1 - h1 * math.tan(theta))
p3 = (p2[0], p2[1] - w2)
p4 = (-p2[0], p2[1])
p5 = (-p3[0], p3[1])

'''Plot'''
figure, ax = plt.subplots()
# ax.set_xlim(left=-11, right=11)
# ax.set_ylim(bottom=-2, top=7)

'''Line endings'''
line1 = [p0, p1]
line2 = [p1, p2]
line3 = [p2, p3]
line4 = [p3, p0]
line5 = [p4, p1]
line6 = [p5, p4]
line7 = [p5, p0]

(line1_xs, line1_ys) = zip(*line1)
(line2_xs, line2_ys) = zip(*line2)
(line3_xs, line3_ys) = zip(*line3)
(line4_xs, line4_ys) = zip(*line4)
(line5_xs, line5_ys) = zip(*line5)
(line6_xs, line6_ys) = zip(*line6)
(line7_xs, line7_ys) = zip(*line7)

'''Create the lines'''
ax.add_line(Line2D(line1_xs, line1_ys, linewidth=1, color='blue'))
ax.add_line(Line2D(line2_xs, line2_ys, linewidth=1, color='red'))
ax.add_line(Line2D(line3_xs, line3_ys, linewidth=1, color='blue'))
ax.add_line(Line2D(line4_xs, line4_ys, linewidth=1, color='red'))
ax.add_line(Line2D(line5_xs, line5_ys, linewidth=1, color='blue'))
ax.add_line(Line2D(line6_xs, line6_ys, linewidth=1, color='red'))
ax.add_line(Line2D(line7_xs, line7_ys, linewidth=1, color='blue'))

'''Show'''
plt.plot()
plt.gca().set_aspect('equal')
plt.show()

'''Horizontal-Tail calculation'''
'''Digital-Model'''
fn = 'GEC-50.igs'
model = Model()

fuselage_width = 2.93
fuselage_height = 3.09
fuselage_len = 11.95

nose_len = 4.66
nose_front_radius = 0.05
dh = 0.6

tail_len = 7.91
tail_up_frame_delta = 0.45
tail_radius = 0.15
tip_back_angle = math.radians(12)

b = fuselage_width / 2
a = fuselage_height / 2

p0 = np.array([0, nose_front_radius, 0])
t0 = np.array([0, 1, 0])
p2 = np.array([nose_len, dh + a, 0])
t2 = np.array([1, 0, 0])
p1c = np.array([nose_len, nose_front_radius, 0])
p1 = p1c + np.array([nose_len * math.cos(math.radians(135)), (p2[1] - nose_front_radius) * math.sin(math.radians(135)), 0])
arc1 = ConicArc(p0, t0, p2, t2, p1)

p3 = np.array([0, -nose_front_radius, 0])
t3 = np.array([0, -1, 0])
p5 = p2 - np.array([0, fuselage_height, 0])
t5 = np.array([1, 0, 0])
p4c = np.array([nose_len, -nose_front_radius, 0])
p4 = p4c + np.array([nose_len * math.cos(math.radians(225)), (-nose_front_radius - p5[1]) * math.sin(math.radians(225)), 0])
arc2 = ConicArc(p3, t3, p5, t5, p4)

p9c = (p2 + p5) / 2
p6 = np.array([0, 0, nose_front_radius])
t6 = (0, 0, 1)
p9 = np.array([0, 0, b]) + p9c
p8 = np.array([p9[0], 0, p9[2]])
t8 = (1, 0, 0)
p7 = np.array([nose_len, 0, nose_front_radius]) + np.array([-nose_len / sqrt2, 0, (b - nose_front_radius) / sqrt2])

arc3 = Arc.from_2pnt(p0, p3, 180, (1, 0, 0))
arc5 = ConicArc(p6, t6, p8, t8, p7)
for i in range(3):
    arc5.Pw[i][1] = i / 2 * p9[1] * arc5.Pw[i][-1]
arc5.reset(arc5.U, arc5.Pw)

'''通过拉伸圆来构造椭圆'''
arc4 = Arc.from_2pnt(p2, p5, 180, (1, 0, 0))
arc4.reset(arc4.U, np.copy(list(map(lambda u: to_homogeneous(to_cartesian(u) * (1, 1, b / a), u[-1]), arc4.Pw))))

arc3_1, arc3_2 = Crv.split(arc3, [0.5])
arc4_1, arc4_2 = Crv.split(arc4, [0.5])

model.add(arc1.to_iges())
model.add(arc2.to_iges())
model.add(arc3_1.to_iges())
model.add(arc3_2.to_iges())
model.add(arc4_1.to_iges())
model.add(arc4_2.to_iges())
model.add(arc5.to_iges())

arc6_1 = deepcopy(arc4_1)
arc6_2 = deepcopy(arc4_2)
arc6_1.pan((fuselage_len, 0, 0))
arc6_2.pan((fuselage_len, 0, 0))
model.add(arc6_1.to_iges())
model.add(arc6_2.to_iges())

line1 = Line(arc4_1(0), arc6_1(0))
line2 = Line(arc4_2(0), arc6_2(0))
line3 = Line(arc4_2(1), arc6_2(1))

model.add(line1.to_iges())
model.add(line2.to_iges())
model.add(line3.to_iges())

fuselage_surf_up = Coons(arc4_1, arc6_1, line1, line2)
model.add(fuselage_surf_up.to_iges())
fuselage_surf_up.mirror('z')
model.add(fuselage_surf_up.to_iges())

fuselage_surf_down = Coons(arc4_2, arc6_2, line2, line3)
model.add(fuselage_surf_down.to_iges())
fuselage_surf_down.mirror('z')
model.add(fuselage_surf_down.to_iges())

p10 = arc6_1(0)
t10 = (1, 0, 0)
p11c = p10 - np.array([0, tail_up_frame_delta, 0])
p12 = p11c + np.array([tail_len, 0, 0])
t12 = (0, -1, 0)
p11 = p11c + np.array([tail_len / sqrt2, tail_up_frame_delta / sqrt2, 0])
arc6 = ConicArc(p10, t10, p12, t12, p11)
model.add(arc6.to_iges())

p13 = arc6_2(1)
t13 = (1, 0, 0)
p15 = p12 - np.array([0, 2 * tail_radius, 0])
t15 = (0, 1, 0)
p14c = np.array([p13[0], p15[1], 0])
p14 = p14c + np.array([tail_len / sqrt2, -(p15[1] - p13[1]) / sqrt2, 0])
arc7 = ConicArc(p13, t13, p15, t15, p14)
model.add(arc7.to_iges())

arc8 = Arc.from_2pnt(p12, p15, 180, (1, 0, 0))
arc8_1, arc8_2 = Crv.split(arc8, [0.5])
model.add(arc8_1.to_iges())
model.add(arc8_2.to_iges())

p16 = arc6_2.start
t16 = (1, 0, 0)
p18 = arc8_2.start
p19 = np.array([p18[0], p16[1], p18[2]])
t19 = (0, 0, -1)
p17c = np.array([p16[0], p16[1], p19[2]])
p17 = p17c + np.array([tail_len / sqrt2, 0, (p16[2] - p19[2]) / sqrt2])

arc9 = ConicArc(p16, t16, p19, t19, p17)
tmp = p18[1] - p16[1]
n = len(arc9.Pw)
for k, pw in enumerate(arc9.Pw):
    w = pw[-1]
    a, b, c = to_cartesian(pw)
    b += k / (n - 1) * tmp
    arc9.Pw[k] = to_homogeneous((a, b, c), w)
arc9.reset(arc9.U, arc9.Pw)
model.add(arc9.to_iges())

'''Wing'''
inner_wing_pan = (7.2, -0.2, 0.8)
foil = ['NACA5312', 'NACA5312']
z_offset = np.array([0., 4.56])
length = np.array([5.8, 5.8])
sweep_back = np.array([0., 0.])
twist = np.array([6., 6.])
dihedral = np.array([0., 0.])
twist_pos = np.array([0.25, 0.25])
y_ref = np.zeros(2)
thickness_factor = np.array([1., 1.])

crv_list = []
for k in range(2):
    wp = WingProfile.from_geom_param(foil[k], z_offset[k], length[k], sweep_back[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness_factor[k])
    crv = wp.crv
    crv.pan(inner_wing_pan)
    crv_list.append(crv)

inner_surf = RuledSurf(crv_list[0], crv_list[1])
model.add(inner_surf.to_iges())
inner_surf.mirror('z')
model.add(inner_surf.to_iges())

inner_tail_surf = Coons(inner_surf.extract('U', 0), inner_surf.extract('U', 1), Line(inner_surf(0, 0), inner_surf(1, 0)), Line(inner_surf(0, 1), inner_surf(1, 1)))
model.add(inner_tail_surf.to_iges())
inner_tail_surf.mirror('z')
model.add(inner_tail_surf.to_iges())

outer_wing_pan = (7.8, -0.1, inner_wing_pan[2] + 4.56)
foil = ['M6', 'M6']
z_offset = np.array([0., 4.48])
length = np.array([3.9, 1.648])
sweep_back = np.array([0., 22.])
twist = np.array([4.5, 4.5])
dihedral = np.array([10., 10.])
twist_pos = np.array([0.25, 0.25])
y_ref = np.zeros(2)
thickness_factor = np.array([1., 1.])

crv_list = []
for k in range(2):
    wp = WingProfile.from_geom_param(foil[k], z_offset[k], length[k], sweep_back[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness_factor[k])
    crv = wp.crv
    crv.pan(outer_wing_pan)
    crv_list.append(crv)

outer_surf = RuledSurf(crv_list[0], crv_list[1])
model.add(outer_surf.to_iges())
outer_surf.mirror('z')
model.add(outer_surf.to_iges())

outer_tail_surf = Coons(outer_surf.extract('U', 0), outer_surf.extract('U', 1), Line(outer_surf(0, 0), outer_surf(1, 0)), Line(outer_surf(0, 1), outer_surf(1, 1)))
model.add(outer_tail_surf.to_iges())
outer_tail_surf.mirror('z')
model.add(outer_tail_surf.to_iges())

'''Vertical tail'''
v_tail_pan = (7.8 + 11.9, 0, 0)
foil = ['NACA0012', 'NACA0012']
z_offset = np.array([0., 2.9988])
length = np.array([3.6812, 1.9151])
sweep_back = np.array([0., 40.])
twist = np.array([0., 0.])
dihedral = np.array([0., 0.])
twist_pos = np.array([0.25, 0.25])
y_ref = np.zeros(2)
thickness_factor = np.array([1., 1.])

crv_list = []
for k in range(2):
    wp = WingProfile.from_geom_param(foil[k], z_offset[k], length[k], sweep_back[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness_factor[k])
    crv = wp.crv
    crv.pan(v_tail_pan)
    crv_list.append(crv)

vtail_surf = RuledSurf(crv_list[0], crv_list[1])
vtail_surf.rotate((0, 0, 0), (-1, 0, 0), 90)
vtail_surf.pan((0, 1.8, 0))
model.add(vtail_surf.to_iges())
vtail_tail_surf = Coons(vtail_surf.extract('U', 0), vtail_surf.extract('U', 1), Line(vtail_surf(0, 0), vtail_surf(1, 0)), Line(vtail_surf(0, 1), vtail_surf(1, 1)))
model.add(vtail_tail_surf.to_iges())

'''Horizontal tail'''
h_tail_pan = (7.8 + 13.3, 4.8, 0)
foil = ['NACA0012', 'NACA0012']
z_offset = np.array([0., 9.742 / 2])
length = np.array([3.2969, 1.7040])
sweep_back = np.array([0., 10.])
twist = np.array([0., 0.])
dihedral = np.array([0., 0.])
twist_pos = np.array([0.25, 0.25])
y_ref = np.zeros(2)
thickness_factor = np.array([1., 1.])

crv_list = []
for k in range(2):
    wp = WingProfile.from_geom_param(foil[k], z_offset[k], length[k], sweep_back[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness_factor[k])
    crv = wp.crv
    crv.pan(h_tail_pan)
    crv_list.append(crv)

htail_surf = RuledSurf(crv_list[0], crv_list[1])
model.add(htail_surf.to_iges())
htail_surf.mirror('z')
model.add(htail_surf.to_iges())

htail_tail_surf = Coons(htail_surf.extract('U', 0), htail_surf.extract('U', 1), Line(htail_surf(0, 0), htail_surf(1, 0)), Line(htail_surf(0, 1), htail_surf(1, 1)))
model.add(htail_tail_surf.to_iges())
htail_tail_surf.mirror('z')
model.add(htail_tail_surf.to_iges())

model.save(fn)
