import unittest
from src.nurbs.spline import Spline
from src.iges.iges_core import IGES_Model
import numpy as np


def build_airfoil(fn):
    fn = '../airfoil/' + fn + '.dat'
    fin = open(fn)
    pnt_list = []
    for pnt in fin:
        x, y, z = pnt.split()
        pnt_list.append([float(x), float(y), 0])

    pts = np.zeros((len(pnt_list), 3))
    for i in range(0, len(pnt_list)):
        for j in range(0, 3):
            pts[i][j] = pnt_list[i][j]

    sp = Spline()
    sp.interpolate(pts)
    return sp


def write_airfoil(airfoil):
    foil = build_airfoil(airfoil)
    model_file = IGES_Model(airfoil + '_spline' + '.igs')
    model_file.add_entity(foil.to_iges())
    model_file.write()


class spline_test(unittest.TestCase):
    def test_airfoil(self):
        write_airfoil('M6')
        write_airfoil('NACA0012')
