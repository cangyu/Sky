import unittest
import numpy as np
from src.nurbs.curve import Arc
from src.iges.iges_core import IGES_Model

try:
    from src.com.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True


class PlanarArcTest(unittest.TestCase):
    radius_list = np.array([0.3, 10], float)
    angle_list = np.array([-20, 0, 3, 45, 65.2, 90, 120, 135, 150, 180, 195, 225, 240, 270, 315, 324, 360, 1024], float)

    @staticmethod
    def build_xy_arc(r, ang):
        crv = Arc(r, ang)
        fn = 'Arc_{}_{}.igs'.format(r, ang)
        arc = IGES_Model(fn)
        arc.add_entity(crv.to_iges(1, 0, [0, 0, 1]))
        arc.write()
        if auto_view:
            view(fn)

    @staticmethod
    def test():
        for radius in PlanarArcTest.radius_list:
            for angle in PlanarArcTest.angle_list:
                PlanarArcTest.build_xy_arc(radius, angle)


class SpatialArcTest(unittest.TestCase):
    @staticmethod
    def test():
        arc1 = Arc(200, 180)
        arc2 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])

        model = IGES_Model('test.igs')
        model.add_entity(arc1.to_iges(1, 0, [0, 0, 1]))
        model.add_entity(arc2.to_iges(1, 0, [0, 1, 0]))
        model.write()

        if auto_view:
            view('test.igs')


if __name__ == '__main__':
    unittest.main()
