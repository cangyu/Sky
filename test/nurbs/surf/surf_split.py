from src.nurbs.surface import ClampedNURBSSurf
from src.iges.iges_core import IGES_Model
import numpy as np

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

U = np.array([0, 0, 0, 0, 1, 1, 1, 1], float)
V = np.array([0, 0, 0, 0.5, 1, 1, 1], float)

alpha_u, beta_u, gamma_u, delta_u = 2, 1, 3, 2
alpha_v, beta_v, gamma_v, delta_v = 2, 1, 3, 2

Pw = np.array([[[0, 0, 0, 1], [0, 1, 1, 1], [0, 2, 3, 1], [0, 3, 2, 1]],
               [[1, 0, 0, 1], [1, 1, 2, 1], [1, 3, 5, 1], [1, 4, 2, 1]],
               [[2, 0, 0, 1], [2, 2, 2, 1], [2, 3, 6, 1], [2, 5, 7, 1]],
               [[3, 0, 0, 1], [3, 1, 1, 1], [3, 2, 2, 1], [3, 3, 3, 1]]], float)

s0 = ClampedNURBSSurf(U, V, Pw)
spt = ClampedNURBSSurf.split(s0, (0.2, 0.3), (0.4, 0.5, 0.6))

fn0 = 'origin.igs'
model0 = IGES_Model(fn0)
model0.add_entity(s0.to_iges())

fn1 = 'split.igs'
model1 = IGES_Model(fn1)
for ue in spt:
    for ve in ue:
        print(ve.U)
        print(ve.V)
        print(ve.Pw)
        model1.add_entity(ve.to_iges())

model0.write()
model1.write()

if auto_view:
    view(fn0)
    view(fn1)
