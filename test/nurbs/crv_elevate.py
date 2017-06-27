import math
import numpy as np
from src.nurbs.curve import NURBS_Curve
from src.iges.iges_core import IGES_Model
from src.com.catia import view

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
C1 = NURBS_Curve(U, Pw)
C1.elevate(1)
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
