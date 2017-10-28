import numpy as np
import math
from copy import deepcopy
from src.aircraft.wing import Wing
from src.aircraft.frame import BWBFrame, chebshev_dist_multi
from src.geom.curve import Line, Arc
from src.geom.surface import BilinearSurf, Coons
from src.msh.tfi import LinearTFI3D
from src.msh.plot3d import PLOT3D, PLOT3D_Block

p3d_msh = PLOT3D()
fn = '3d_multi_blk_grid.xyz'

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

C10 = deepcopy(C01)
C10.reverse()
C98 = deepcopy(C89)
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

C32 = Arc.from_2pnt(P[3], P[2], 180, [0, 0, -1])
C1110 = Arc.from_2pnt(P[11], P[10], 180, [0, 0, -1])
C23 = Arc.from_2pnt(P[2], P[3], 180, [0, 0, 1])
C1011 = Arc.from_2pnt(P[10], P[11], 180, [0, 0, 1])

S1 = Coons(C32, C1110, Line(P[3], P[11]), Line(P[2], P[10]))
S3 = Coons(Line(P[3], P[11]), C19, Line(P[3], P[1]), Line(P[11], P[9]))
S4 = Coons(Line(P[2], P[10]), C08, Line(P[2], P[0]), Line(P[10], P[8]))
S5 = Coons(Line(P[3], P[1]), Line(P[2], P[0]), C32, C10)
S6 = Coons(Line(P[11], P[9]), Line(P[10], P[8]), C1110, C98)

blk0_tfi_grid = LinearTFI3D(lambda v, w: S1(v, w),
                            lambda v, w: sf(1.0 - v, w),
                            lambda w, u: S3(w, u),
                            lambda w, u: S4(w, u),
                            lambda u, v: S5(u, v),
                            lambda u, v: S6(u, v))

U, V, W = 32, 42, 25
u_list = np.linspace(0, 1.0, U + 1)
v_list = np.linspace(0, 1.0, V + 1)
w_list = np.linspace(0, 1.0, W + 1)
blk0_tfi_grid.calc_grid(u_list, v_list, w_list)
p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk0_tfi_grid.get_grid()))

SS1 = Coons(Line(P[0], P[2]), Line(P[8], P[10]), C08, Line(P[2], P[10]))
SS2 = BilinearSurf(np.array([[P[4], P[12]], [P[6], P[14]]]))
SS3 = Coons(C08, Line(P[4], P[12]), Line(P[0], P[4]), Line(P[8], P[12]))
SS4 = BilinearSurf(np.array([[P[2], P[6]], [P[10], P[14]]]))
SS5 = BilinearSurf(np.array([[P[0], P[2]], [P[4], P[6]]]))
SS6 = BilinearSurf(np.array([[P[8], P[10]], [P[12], P[14]]]))

blk1_tfi_grid = LinearTFI3D(lambda v, w: SS1(v, w),
                            lambda v, w: SS2(v, w),
                            lambda w, u: SS3(w, u),
                            lambda w, u: SS4(w, u),
                            lambda u, v: SS5(u, v),
                            lambda u, v: SS6(u, v))

N = 15
n_list = np.linspace(0.0, 1.0, N + 1)
blk1_tfi_grid.calc_grid(n_list, u_list, w_list)
p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk1_tfi_grid.get_grid()))

SSS1 = Coons(Line(P[3], P[1]), Line(P[11], P[9]), Line(P[3], P[11]), C19)
SSS2 = BilinearSurf(np.array([[P[7], P[15]], [P[5], P[13]]]))
SSS3 = BilinearSurf(np.array([[P[3], P[7]], [P[11], P[15]]]))
SSS4 = Coons(C19, Line(P[5], P[13]), Line(P[1], P[5]), Line(P[9], P[13]))
SSS5 = BilinearSurf(np.array([[P[3], P[1]], [P[7], P[5]]]))
SSS6 = BilinearSurf(np.array([[P[11], P[9]], [P[15], P[13]]]))

blk2_tfi_grid = LinearTFI3D(lambda v, w: SSS1(v, w),
                            lambda v, w: SSS2(v, w),
                            lambda w, u: SSS3(w, u),
                            lambda w, u: SSS4(w, u),
                            lambda u, v: SSS5(u, v),
                            lambda u, v: SSS6(u, v))

blk2_tfi_grid.calc_grid(n_list, u_list, w_list)
p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk2_tfi_grid.get_grid()))

SSSI1 = Coons(Line(P[0], P[1]), Line(P[8], P[9]), C08, C19)
SSSI2 = BilinearSurf(np.array([[P[4], P[12]], [P[5], P[13]]]))
SSSI3 = Coons(C08, Line(P[4], P[12]), Line(P[0], P[4]), Line(P[8], P[12]))
SSSI4 = Coons(C19, Line(P[5], P[13]), Line(P[1], P[5]), Line(P[9], P[13]))
SSSI5 = BilinearSurf(np.array([[P[0], P[1]], [P[4], P[5]]]))
SSSI6 = BilinearSurf(np.array([[P[8], P[9]], [P[12], P[13]]]))

blk3_tfi_grid = LinearTFI3D(lambda v, w: SSSI1(v, w),
                            lambda v, w: SSSI2(v, w),
                            lambda w, u: SSSI3(w, u),
                            lambda w, u: SSSI4(w, u),
                            lambda u, v: SSSI5(u, v),
                            lambda u, v: SSSI6(u, v))

M = 6
m_list = np.linspace(0, 1, M + 1)
blk3_tfi_grid.calc_grid(n_list, m_list, w_list)
p3d_msh.add_block(PLOT3D_Block.build_from_3d(blk3_tfi_grid.get_grid()))
p3d_msh.write(fn)
