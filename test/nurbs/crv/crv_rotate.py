import unittest
from src.nurbs.curve import Arc
from src.iges.iges_core import IGES_Model

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True


class RotateTest(unittest.TestCase):
    @staticmethod
    def test():
        arc1 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
        arc2 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
        arc2.rotate([0, 0, 0], [1, 1, 1], 30)
        arc3 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
        arc3.rotate([0, 0, 0], [0, 0, 1], 45)
        arc4 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
        arc4.rotate([0, 0, 0], [0, 0, 1], 90)
        arc5 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
        arc5.rotate([0, 0, 0], [0, 0, 1], 120)

        model = IGES_Model('test.igs')
        model.add_entity(arc1.to_iges())
        model.add_entity(arc2.to_iges())
        model.add_entity(arc3.to_iges())
        model.add_entity(arc4.to_iges())
        model.add_entity(arc5.to_iges())
        model.write()

        if auto_view:
            view('test.igs')


if __name__ == '__main__':
    unittest.main()
