import math
import numpy as np
from src.nurbs.curve import NURBS_Curve
from src.iges.iges_core import IGES_Model
from src.com.catia import view

auto_view = False

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
C0 = NURBS_Curve(U, Pw)
model = IGES_Model('before.igs')
model.add_entity(C0.to_iges())
model.write()

print('C0:')
print('Knot vector:')
print(C0.U)
print('Control points:')
print(C0.Pw)

'''插入1个不在原节点矢量中的节点'''
C1 = NURBS_Curve(U, Pw)
C1.insert_knot(2.5)
model = IGES_Model('after1.igs')
model.add_entity(C1.to_iges())
model.write()

print('C1:')
print('Knot vector:')
print(C1.U)
print('Control points:')
print(C1.Pw)

'''插入1个在原节点矢量中的节点'''
C2 = NURBS_Curve(U, Pw)
C2.insert_knot(2)
model = IGES_Model('after2.igs')
model.add_entity(C2.to_iges())
model.write()

print('C2:')
print('Knot vector:')
print(C2.U)
print('Control points:')
print(C2.Pw)

'''插入1个不在原节点矢量中的节点3次'''
C3 = NURBS_Curve(U, Pw)
C3.insert_knot(2.5, 3)
model = IGES_Model('after3.igs')
model.add_entity(C3.to_iges())
model.write()

print('C3:')
print('Knot vector:')
print(C3.U)
print('Control points:')
print(C3.Pw)

'''插入1个在原节点矢量中的节点2次'''
C4 = NURBS_Curve(U, Pw)
C4.insert_knot(2, 2)
model = IGES_Model('after4.igs')
model.add_entity(C4.to_iges())
model.write()

print('C4:')
print('Knot vector:')
print(C4.U)
print('Control points:')
print(C4.Pw)

'''逐次插入1个在原节点矢量中的节点2次'''
C5 = NURBS_Curve(U, Pw)
C5.insert_knot(2)
C5.insert_knot(2)
model = IGES_Model('after5.igs')
model.add_entity(C5.to_iges())
model.write()

print('C5:')
print('Knot vector:')
print(C5.U)
print('Control points:')
print(C5.Pw)

'''插入1个不在原节点矢量中的节点4次'''
C6 = NURBS_Curve(U, Pw)
try:
    C6.insert_knot(2.5, 4)
except ValueError as ve:
    print('Ok, illegal insertion detected with following msg:')
    print(ve)
else:
    print('WTF?')

if auto_view:
    view('before.igs')
    view('after1.igs')
    view('after2.igs')
    view('after3.igs')
    view('after4.igs')
    view('after5.igs')
