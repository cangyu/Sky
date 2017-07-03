import math
import numpy as np
from src.nurbs.curve import NURBS_Curve
from src.iges.iges_core import IGES_Model
from src.com.catia import view

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

X = [2.5, 2.5, 2.5]

'''原始曲线'''
C0 = NURBS_Curve(U, Pw)
model = IGES_Model('before.igs')
model.add_entity(C0.to_iges())
model.write()

'''插入若干节点'''
C1 = NURBS_Curve(U, Pw)
C1.refine(X)
model = IGES_Model('after.igs')
model.add_entity(C1.to_iges())
model.write()

print('Original Knot vector:')
print(C0.U)
print('Original Control points:')
print(C0.Pw)

print('New Knot vector:')
print(C1.U)
print('New Control points:')
print(C1.Pw)

if auto_view:
    view('before.igs')
    view('after.igs')
