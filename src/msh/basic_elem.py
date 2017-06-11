import numpy as np


class Point(object):
    def __init__(self, index, *args):
        self.index = index
        self.dim = len(args)
        self.x = np.empty(self.dim, float)
        for i in range(self.dim):
            self.x[i] = args[i]


class Face(object):
    def __init__(self, index, *args):
        self.index = index
        self.pts = []
        for pnt in args:
            self.pts.append(pnt)


class Cell(object):
    def __init__(self, index, *args):
        self.index = index
        self.face_list = []
        for face in args:
            self.face_list.append(face)


class Mesh(object):
    def __init__(self):
        self.pnt_list = []
        self.face_list = []
        self.cell_list = []
