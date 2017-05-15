import math
import os
import numpy as np
from src.nurbs.curve import *
from src.iges.iges_entity126 import *
from src.iges.iges_entity110 import *
from src.iges.iges_entity116 import *
from src.iges.iges_core import *

AIRFOIL_DIR = '../airfoil/'
AIRFOIL_LIST = []


def update_airfoil_list():
    for f in os.listdir(AIRFOIL_DIR):
        base, ext = os.path.splitext(f)
        if ext == '.dat':
            AIRFOIL_LIST.append(base)


class Airfoil(object):
    '''
    2D Airfoil, chord length equals to 1
    '''

    def __init__(self, _filename):
        self.name = _filename
        self.x = []
        self.y_up = []
        self.y_down = []

        # Read input data
        airfoil = open(AIRFOIL_DIR + self.name + ".dat")
        for line in airfoil:
            (_x, _y_up, _y_down) = line.split()
            self.x.append(float(_x))
            self.y_up.append(float(_y_up))
            self.y_down.append(float(_y_down))
        airfoil.close()

        self.n = len(self.x)
        self.pts_num = 2 * self.n - 1
        self.pts = np.zeros((self.pts_num, 3), float)

        # Arrange in XFoil's format
        cnt = 0
        for i in range(0, self.n):
            self.pts[cnt][0] = self.x[- 1 - i]
            self.pts[cnt][1] = self.y_up[- 1 - i]
            cnt += 1

        for i in range(1, self.n):
            self.pts[cnt][0] = self.x[i]
            self.pts[cnt][1] = self.y_down[i]
            cnt += 1

    def get_pts(self):
        return self.pts

    def iges(self):
        foil = IGES_Model(self.name + '.igs')

        # Points
        for i in range(0, 2 * self.n - 1):
            foil.AddPart(IGES_Entity116(self.pts[i][0], self.pts[i][1], self.pts[i][2]))

        # Closed Curve
        closed_pts = np.zeros((self.pts_num + 1, 3), float)
        for j in range(0, 3):
            closed_pts[self.pts_num][j] = self.pts[0][j]
        for i in range(0, self.pts_num):
            for j in range(0, 3):
                closed_pts[i][j] = self.pts[i][j]
        foil.AddPart(Spline(closed_pts).iges())

        foil.Generate()


class Wing_Profile(object):
    '''
    3D profile for a physical wing at given position
    '''

    def __init__(self, _airfoil, _ends, _thickness_factor=1.0):

        self.airfoil = Airfoil(_airfoil)
        self.ends = _ends
        self.thickness = _thickness_factor
        self.n = self.airfoil.n

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
            self.pts[0][0][i] = self.pts[1][0][i] = float(chord_len * self.airfoil.x[i])

            # y-dir
            self.pts[0][1][i] = float(chord_len * self.airfoil.y_up[i])
            self.pts[1][1][i] = float(chord_len * self.airfoil.y_down[i])

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

    def getPointList(self):
        ret = np.zeros((2 * self.n - 1, 3), float)
        cpi = 0

        for i in range(0, self.n):
            ret[cpi][0] = self.pts[0][0][self.n - 1 - i]
            ret[cpi][1] = self.pts[0][1][self.n - 1 - i]
            ret[cpi][2] = self.pts[0][2][self.n - 1 - i]
            cpi += 1

        for i in range(1, self.n):
            ret[cpi][0] = self.pts[1][0][i]
            ret[cpi][1] = self.pts[1][1][i]
            ret[cpi][2] = self.pts[1][2][i]
            cpi += 1

        return ret
