import copy
import math

import numpy as np
from src.msh.plot3d import PLOT3D, PLOT3D_Block
from src.msh.tfi import LinearTFI2D, LinearTFI3D

from nurbs import Line, Circle
from src.aircraft.frame import BWBFrame, chebshev_dist_multi
from src.aircraft.wing import Wing

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True


def write_uniform_p3d(msh, u, v, w, _fn="msh_p3d.xyz"):
    """
    构建单块Plot3D网格
    :param msh: 3D TFI Grid
    :type msh: LinearTFI3D
    :param u: U方向网格数量
    :type u: int
    :param v: V方向网格数量
    :type v: int
    :param w: W方向网格数量
    :type w: int
    :param _fn: 输出文件名
    :type _fn: str
    :return: None
    """

    u_list = np.linspace(0, 1.0, u + 1)
    v_list = np.linspace(0, 1.0, v + 1)
    w_list = np.linspace(0, 1.0, w + 1)
    msh.calc_grid(u_list, v_list, w_list)
    grid = PLOT3D()
    grid.add_block(PLOT3D_Block.build_from_3d(msh.get_grid()))
    grid.write(_fn)


c_root = 10
c_mid = 6
c_tip = 1.1
b_mid = 4.5
b_tip = 13
alpha_mid = 35
alpha_tip = 30
frm = BWBFrame(c_root, c_mid, c_tip, b_mid, b_tip, alpha_mid, alpha_tip)

inner_sec_num = 6
outer_sec_num = 8
u_mid = b_mid / b_tip
u_dist = chebshev_dist_multi([0, u_mid, 1], [inner_sec_num, outer_sec_num])
sec_num = len(u_dist)

foil = ['NACA0012', 'NACA0012', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6']
z_offset = list(map(frm.z, u_dist))
length = list(map(lambda u: frm.x_tail(u) - frm.x_front(u), u_dist))
sweep_back = list(map(lambda u: math.degrees(math.atan2(frm.x_front(u), frm.z(u))), u_dist))
twist = np.array([3, 2.5, 2.2, 2, 2, 1.8, 1.6, 1.1, 1, 1, 0.5, 0.2, 0])
dihedral = np.array([0, 1, 1.1, 1.2, 2, 2, 2, 2, 2, 2, 2.5, 3, 3])
twist_pos = np.full(sec_num, 0.25)
y_ref = np.zeros(sec_num)
thickness_factor = np.ones(sec_num)

wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)

La = length[0]
Lt = 40 * La
R = 20 * La
Span = z_offset[-1]

sf = wg.surf()
C01 = Line(wg.section[0].tail_up, wg.section[0].tail_down)
C89 = Line(wg.section[-1].tail_up, wg.section[-1].tail_down)
C08 = sf.extract('U', 0)
C19 = sf.extract('U', 1)

C10 = copy.deepcopy(C01)
C10.reverse()
C98 = copy.deepcopy(C89)
C98.reverse()

P = np.zeros((16, 3))
P[0] = sf(0, 0)
P[1] = sf(1, 0)
P[8] = sf(0, 1)
P[9] = sf(1, 1)
P[2] = np.copy(P[0])
P[2][1] += R
P[3] = np.copy(P[1])
P[3][1] -= R
P[10] = np.copy(P[8])
P[10][1] += R
P[11] = np.copy(P[9])
P[11][1] -= R
P[4] = np.array([La + Lt, P[0][1], 0.0], float)
P[5] = np.array([La + Lt, P[1][1], 0.0], float)
P[6] = np.copy(P[4])
P[6][1] += R
P[7] = np.copy(P[5])
P[7][1] -= R
P[12] = np.copy(P[8])
P[12][0] = P[4][0]
P[13] = np.copy(P[9])
P[13][0] = P[5][0]
P[14] = np.copy(P[12])
P[14][1] = P[6][1]
P[15] = np.copy(P[13])
P[15][1] = P[7][1]

LineEnding = np.array([[2, 0], [1, 0], [3, 1], [10, 8], [9, 8], [11, 9], [3, 11], [2, 10], [2, 6], [0, 4],
                       [1, 5], [3, 7], [10, 14], [8, 12], [9, 13], [11, 15], [4, 6], [5, 4], [7, 5], [12, 14], [13, 12], [15, 13], [6, 14], [4, 12], [5, 13], [7, 15]])
LineName = []
for i in range(len(LineEnding)):
    LineName.append("L{}{}".format(LineEnding[i][0], LineEnding[i][1]))

LineMap = {}
for i in range(len(LineName)):
    LineMap[LineName[i]] = i

L = []
for ln in LineName:
    sp, ep = LineEnding[LineMap[ln]]
    L.append(Line(P[sp], P[ep]))

C32 = Circle.from_2pnt(P[3], P[2], 180, [0, 0, -1])
C1110 = Circle.from_2pnt(P[11], P[10], 180, [0, 0, -1])
C23 = Circle.from_2pnt(P[2], P[3], 180, [0, 0, 1])
C1011 = Circle.from_2pnt(P[10], P[11], 180, [0, 0, 1])

'''
fn = "frame.igs"
frame = IGES_Model(fn)
for line in L:
    frame.add_entity(line.to_iges())
frame.add_entity(C08.to_iges(0, 0, [0, 0, 0]))
frame.add_entity(C19.to_iges(0, 0, [0, 0, 0]))
frame.add_entity(C01.to_iges())
frame.add_entity(C89.to_iges())
frame.add_entity(C23.to_iges(1, 0, [0, 0, 1]))
frame.add_entity(C1011.to_iges(1, 0, [0, 0, 1]))
frame.write()
if auto_view:
    view(fn)
'''

S1 = LinearTFI2D(C32, L[LineMap['L311']], C1110, L[LineMap['L210']])
S3 = LinearTFI2D(L[LineMap['L311']], L[LineMap['L31']], C19, L[LineMap['L119']])
S4 = LinearTFI2D(L[LineMap['L210']], L[LineMap['L20']], C08, L[LineMap['L108']])
S5 = LinearTFI2D(L[LineMap['L31']], C32, L[LineMap['L20']], C10)
S6 = LinearTFI2D(L[LineMap['L119']], C1110, L[LineMap['L108']], C98)

tfi_grid = LinearTFI3D(lambda v, w: S1(v, w),
                       lambda v, w: sf(1.0 - v, w),
                       lambda w, u: S3(w, u),
                       lambda w, u: S4(w, u),
                       lambda u, v: S6(u, v),
                       lambda u, v: S6(u, v))

write_uniform_p3d(tfi_grid, 45, 28, 16)
