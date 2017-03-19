import sys
import os
import math
import time
import numpy as np
import pylab as pl
from abc import ABCMeta, abstractmethod
from src.iges.iges_core import IGES_Model
from src.iges.iges_entity110 import *
from src.iges.iges_entity112 import *
from src.iges.iges_entity116 import *


class Airfoil(object):
    '''
    2D Airfoil, with chord length equals 1
    '''

    def __init__(self, _filename):
        self.x = []
        self.y_up = []
        self.y_down = []

        airfoil = open(_filename)
        for line in airfoil:
            (_x, _y_up, _y_down) = line.split()
            self.x.append(float(_x))
            self.y_up.append(float(_y_up))
            self.y_down.append(float(_y_down))
        airfoil.close()


class Wing_Profile(object):
    '''
    3D profile for a physical wing at given position
    '''

    def __init__(self, _airfoil, _ends, _thickness_factor=1.0):

        self.airfoil = _airfoil
        self.ends = _ends
        self.thickness = _thickness_factor
        self.n = len(_airfoil.x)

        assert _ends[0][2] == _ends[1][2]

        chord_len = math.sqrt(math.pow(_ends[0][0] - _ends[1][0], 2) + math.pow(_ends[0][1] - _ends[1][1], 2))
        assert chord_len > 0

        dx = _ends[1][0] - _ends[0][0]
        dy = _ends[1][1] - _ends[0][1]
        dx /= chord_len
        dy /= chord_len
        rotation = complex(dx, dy)

        # 2 line, 3 dimension, n point on each line
        self.pts = np.zeros((2, 3, self.n), dtype=float)

        # stretch
        for i in range(0, self.n):
            # x-dir
            self.pts[0][0][i] = self.pts[1][0][i] = float(chord_len * _airfoil.x[i])

            # y-dir
            self.pts[0][1][i] = float(chord_len * _airfoil.y_up[i])
            self.pts[1][1][i] = float(chord_len * _airfoil.y_down[i])

            # z-dir
            self.pts[0][2][i] = self.pts[1][2][i] = _ends[0][2]

        # thickness
        for i in range(0, self.n):
            self.pts[0][1][i] *= _thickness_factor
            self.pts[1][1][i] *= _thickness_factor

        # rotate
        for k in range(0, 2):
            for i in range(0, self.n):
                ori_vect = complex(self.pts[k][0][i], self.pts[k][1][i])
                ori_vect *= rotation

                self.pts[k][0][i] = ori_vect.real
                self.pts[k][1][i] = ori_vect.imag

        # move to target
        for k in range(0, 2):
            for i in range(0, self.n):
                self.pts[k][0][i] += _ends[0][0]
                self.pts[k][1][i] += _ends[0][1]

    def AttachTo(self, _model):
        _model.AddPart(IGES_Entity112_Builder(self.pts[0][0], self.pts[0][0], self.pts[0][1], self.pts[0][2], ).GetEntity())
        _model.AddPart(IGES_Entity112_Builder(self.pts[1][0], self.pts[1][0], self.pts[1][1], self.pts[1][2], ).GetEntity())


z = np.array([0., 0.624029, 1.38967, 2.43503, 3.73439, 5.25574, 6.96162,
              8.81003, 10.7555, 12.75, 14.7445, 16.69, 18.5384, 20.2443, 21.7656,
              23.065, 24.1103, 24.876, 25.343, 25.5])

x_front = np.array([0., 0.05, 0.3, 1.7, 4.9, 6.85, 8.45, 9.65, 10.6, 11.1, 11.7, 12.1,
                    12.4, 12.8, 13.2, 13.7, 14.1, 14.5, 15.2, 16.])

y_front = np.linspace(0, 0, len(z))

x_tail = np.array([19.7, 19.6, 19.6, 19.5, 19.3, 19, 18.3, 17.3, 16.6, 16.5,
                   16.8, 17, 17.45, 17.8, 18.1, 18.4, 18.55, 18.65, 18.3, 17.8])

y_tail = np.linspace(0, 0, len(z))

root_line = np.array([[x_front[0], y_front[0], z[0]],
                      [x_tail[0], y_tail[0], z[0]]])

tip_line = np.array([[x_front[len(z) - 1], y_front[len(z) - 1], z[len(z) - 1]],
                     [x_tail[len(z) - 1], y_tail[len(z) - 1], z[len(z) - 1]]])

plane = IGES_Model()

plane.AddPart(IGES_Entity112_Builder(z, x_front, y_front, z).GetEntity())
plane.AddPart(IGES_Entity112_Builder(z, x_tail, y_tail, z).GetEntity())

plane.AddPart(IGES_Entity110_Builder(root_line).GetEntity())
plane.AddPart(IGES_Entity110_Builder(tip_line).GetEntity())

for i in range(0, len(z)):
    plane.AddPart(IGES_Entity116_Builder(x_front[i], y_front[i], z[i]).GetEntity())
    plane.AddPart(IGES_Entity116_Builder(x_tail[i], y_tail[i], z[i]).GetEntity())

naca0012 = Airfoil("../airfoil/naca0012.dat")


wp0=Wing_Profile(naca0012, root_line)
wp1=Wing_Profile(naca0012, tip_line)

wp0.AttachTo(plane)
wp1.AttachTo(plane)

plane.Generate()
