import math
import numpy as np
from src.nurbs.utility import equal
from src.nurbs.curve import NURBS_Curve
from src.iges.iges_core import IGES_Model

try:
    from src.com.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = False

U = np.array([0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 1, 1, 1, 1], float)
Pw = np.array([[0, 0, 0, 1],
               [0, math.pi, 0, 1],
               [0, 0, 4, 1],
               [1, 0, -2, 1],
               [0, 1, 0, 1],
               [2, 0, 0, 1],
               [0, 0, 9, 1],
               [0.618, 1.414, 2.718, 1]], float)

'''原始曲线'''
C0 = NURBS_Curve(U, Pw)
model = IGES_Model('before.igs')
model.add_entity(C0.to_iges())
model.write()

print('Original Knot vector:')
print(C0.U)
print('Original Control points:')
print(C0.Pw)

'''提升2阶'''
C1 = C0.elevate(2, return_raw=False)
model = IGES_Model('after.igs')
model.add_entity(C1.to_iges())
model.write()

print('New Knot vector:')
print(C1.U)
print('New Control points:')
print(C1.Pw)

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
