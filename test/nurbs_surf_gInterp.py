import unittest
import numpy as np
from src.iges.iges_core import IGES_Model
from src.nurbs.nurbs_surface import GlobalInterpolatedSurf


def read_airfoil_pts(fn):
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

    return pts


def write_wing(fn, N, L, p, q):
    N += 1
    pts = read_airfoil_pts(fn)
    m, dim = pts.shape
    all_pts = np.zeros((N, m, dim))
    for i in range(0, N):
        all_pts[i] = np.copy(pts)
        for j in range(0, m):
            all_pts[i][j][-1] = L * i

    wsf = GlobalInterpolatedSurf(all_pts, p, q)
    model_file = IGES_Model("{}_{}_{}_{}.igs".format(fn, p, q, (N - 1) * L))
    model_file.add_entity(wsf.to_iges(0, 0, 0, 0))
    model_file.write()


class nurbs_surf_gInterp_test(unittest.TestCase):
    def test_wing(self):
        write_wing('M6', 10, 0.7, 3, 5)
        write_wing('M6', 10, 0.7, 5, 5)
        write_wing('NACA0012', 10, 0.7, 3, 5)
        write_wing('NACA0012', 10, 0.7, 5, 5)
