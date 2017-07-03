import unittest
from src.nurbs.curve import Arc
from src.nurbs.surface import ExtrudedSurf
from src.iges.iges_core import IGES_Model

try:
    from src.com.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True


def build_extrued(crv, direction):
    """
    根据给定曲线与方向构造拉伸曲面
    :param crv: 基准曲线
    :param direction: 拉伸方向矢量
    :return: None.
    """

    surf = ExtrudedSurf(crv, direction)
    fn = 'test.igs'
    model_file = IGES_Model(fn)
    model_file.add_entity(surf.to_iges(0, 0, 0, 0))
    model_file.write()

    if auto_view:
        view(fn)


class ExtrudedSurfTest(unittest.TestCase):
    @staticmethod
    def test():
        build_extrued(Arc.from_2pnt([30, 45, 69], [44, 66, 88], 75, [1, -1, 1]), [30, 30, 30])
        build_extrued(Arc.from_2pnt([0, 50, 0], [50, 0, 0], 315, [0, 0, 1]), [0, 0, 90])


if __name__ == '__main__':
    unittest.main()
