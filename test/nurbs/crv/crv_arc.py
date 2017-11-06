import unittest
import numpy as np
from nurbs import Arc
from iges import Model


class PlanarArcTest(unittest.TestCase):
    radius_list = np.array([0.3, 10], float)
    angle_list = np.array([-20, 0, 3, 45, 65.2, 90, 120, 135, 150, 180, 195, 225, 240, 270, 315, 324, 360, 1024], float)

    @staticmethod
    def build_xy_arc(r, ang):
        crv = Arc(r, ang)
        fn = 'Arc_{}_{}.igs'.format(r, ang)
        arc = Model()
        arc.add(crv.to_iges(1, 0, [0, 0, 1]))
        arc.save(fn)

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

        model = Model()
        model.add(arc1.to_iges(1, 0, [0, 0, 1]))
        model.add(arc2.to_iges(1, 0, [0, 1, 0]))
        model.save('test.igs')


if __name__ == '__main__':
    unittest.main()
