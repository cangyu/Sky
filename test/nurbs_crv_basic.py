import unittest
import numpy as np
import math
from src.nurbs.utility import to_homogeneous
from src.nurbs.curve import NURBS_Curve
from src.iges.iges_core import IGES_Model


def build_nurbs_crv_2D(U, w, P, z=0):
    Pw = np.zeros((len(P), 4))
    for i in range(0, len(P)):
        temp = np.array([P[i][0], P[i][1], z])
        Pw[i] = to_homogeneous(temp, w[i])

    return NURBS_Curve(U, Pw)


class nurbs_curve_test(unittest.TestCase):
    def test_circle1(self):
        U = np.array([0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1])
        w = np.array([1, 1 / math.sqrt(2), 1, 1 / math.sqrt(2), 1, 1 / math.sqrt(2), 1, 1 / math.sqrt(2), 1])
        P = np.array([[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1], [1, 0]])

        circle = build_nurbs_crv_2D(U, w, P)
        iges_file = IGES_Model('circle1.igs')
        iges_file.add_entity(circle.to_iges(1, 0, [0, 0, 1]))
        iges_file.write()

    def test_circle2(self):
        U = np.array([0, 0, 0, 0.5, 1, 1, 1])
        w = np.array([1, 0.5, 0.5, 1])
        P = np.array([[1, 0], [1, 1], [-1, 1], [-1, 0]])

        circle = build_nurbs_crv_2D(U, w, P)
        iges_file = IGES_Model('circle2.igs')
        iges_file.add_entity(circle.to_iges(1, 0, [0, 0, 1]))
        iges_file.write()

    def test_circle3(self):
        U = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        w = np.array([1, 1 / 3, 1 / 3, 1])
        P = np.array([[1, 0], [1, 2], [-1, 2], [-1, 0]])

        circle = build_nurbs_crv_2D(U, w, P)
        iges_file = IGES_Model('circle3.igs')
        iges_file.add_entity(circle.to_iges(1, 0, [0, 0, 1]))
        iges_file.write()

    def test_circle4(self):
        a = math.sqrt(3) / 2
        U = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        w = np.array([1, 1 / 6, 1 / 6, 1])
        P = np.array([[a, 1 / 2], [2 * a, -3], [-2 * a, -3], [-a, 1 / 2]])

        circle = build_nurbs_crv_2D(U, w, P)
        iges_file = IGES_Model('circle4.igs')
        iges_file.add_entity(circle.to_iges(1, 0, [0, 0, 1]))
        iges_file.write()

    def test_line(self):
        U = np.array([0, 0, 1, 1])
        w = np.array([1, 1])
        P = np.array([[10, 10], [100, 100]])

        line = build_nurbs_crv_2D(U, w, P)
        iges_file = IGES_Model('line.igs')
        iges_file.add_entity(line.to_iges(1, 0, [0, 0, 1]))
        iges_file.write()
