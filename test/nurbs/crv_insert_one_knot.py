import math
import numpy as np
from src.nurbs.curve import NURBS_Curve
from src.iges.iges_core import IGES_Model
from src.com.catia import view

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

'''插入1个不在原节点矢量中的节点'''
C1 = NURBS_Curve(U, Pw)
C1.insert_one_knot(2.5)
model = IGES_Model('after1.igs')
model.add_entity(C1.to_iges())
model.write()

'''插入1个在原节点矢量中的节点'''
C2 = NURBS_Curve(U, Pw)
C2.insert_one_knot(2)
model = IGES_Model('after2.igs')
model.add_entity(C2.to_iges())
model.write()

view('before.igs')
view('after1.igs')
view('after2.igs')
