import math
import numpy as np
from src.nurbs.curve import ClampedNURBSCrv
from src.nurbs.utility import equal
from src.iges.iges_core import IGES_Model

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

U = np.array([0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5], float)
Pw = np.array([[0, 0, 0, 1],
               [0, math.pi, 0, 1],
               [0, 0, 4, 1],
               [1, 0, -2, 1],
               [0, 1, 0, 1],
               [2, 0, 0, 1],
               [0, 0, 9, 1],
               [0.618, 1.414, 2.718, 1]], float)

'''原始曲线'''
C0 = ClampedNURBSCrv(U, Pw)
model = IGES_Model('before.igs')
model.add_entity(C0.to_iges())
model.write()

print('C0:')
print('Knot vector:')
print(C0.U)
print('Control points:')
print(C0.Pw)

'''Bezier Decompose'''
C1 = ClampedNURBSCrv.decompose(C0)
model = IGES_Model('after.igs')
for crv in C1:
    model.add_entity(crv.to_iges())
model.write()

if auto_view:
    view('before.igs')
    view('after.igs')

wtf = False
N = 1000
u = np.linspace(0.0, 1.0, N)
v0 = np.empty((N, 3), float)
v1 = np.empty((N, 3), float)
for i in range(N):
    v0[i] = np.copy(C0(u[i]))
    v1[i] = np.copy(C1(u[i]))
    tmp = v0[i] - v1[i]
    if not equal(np.linalg.norm(tmp, 2), 0.0):
        print(i)
        print(v0[i])
        print(v1[i])
        if not wtf:
            wtf = True

if wtf:
    print('WTF?')
else:
    print('OK, all close!')
