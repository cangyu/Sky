import unittest
import numpy as np
import math
from src.nurbs.curve import Arc, NURBS_Curve
from src.iges.iges_core import IGES_Model

angle_list = np.array([-20, 0, 3, 45, 65.2, 90, 120, 135, 150, 180, 195, 225, 240, 270, 315, 324, 360, 1024], float)
radius_list = np.array([0.3, 10], float)


def build_xy_arc(r, ang):
    U, Pw = Arc.build_simplified_arc(r, ang)
    crv = NURBS_Curve(U, Pw)
    model = IGES_Model('Arc_{}_{}.igs'.format(r, ang))
    model.add_entity(crv.to_iges(1, 0, [0, 0, 1]))
    model.write()


class SimpleConicArcTest(unittest.TestCase):
    @staticmethod
    def test_arc():
        for radius in radius_list:
            for angle in angle_list:
                build_xy_arc(radius, angle)
