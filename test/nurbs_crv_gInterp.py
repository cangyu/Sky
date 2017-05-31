import unittest
import numpy as np
from src.nurbs.curve import GlobalInterpolatedCrv
from src.iges.iges_core import IGES_Model


def build_airfoil(fn, p=3):
    fn = '../airfoil/' + fn + '.dat'
    pnt_list = []
    fin = open(fn)
    for pnt in fin:
        x, y, z = pnt.split()
        pnt_list.append([float(x), float(y), 0])
    fin.close()

    pts = np.zeros((len(pnt_list), 3))
    for i in range(0, len(pnt_list)):
        for j in range(0, 3):
            pts[i][j] = pnt_list[i][j]

    return GlobalInterpolatedCrv(pts, p)


def write_airfoil(airfoil, p):
    foil = build_airfoil(airfoil, p)
    model_file = IGES_Model(airfoil + '_' + str(p) + '.igs')
    model_file.add_entity(foil.to_iges(1, 0, [0, 0, 1]))
    model_file.write()


class nurbs_curve_test(unittest.TestCase):
    def test_airfoil(self):
        write_airfoil('M6', 3)
        write_airfoil('M6', 5)
        write_airfoil('NACA0012', 3)
        write_airfoil('NACA0012', 5)
