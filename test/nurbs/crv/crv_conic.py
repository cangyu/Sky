import math

from nurbs import ConicArc
from src.iges import Model

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True


def show(crv, fn):
    """
    将曲线输出并显示
    :param crv: Curve to be inspected
    :type crv: Crv
    :param fn: File name
    :type fn: str
    :return: None
    """

    _model = Model()
    _model.add_entity(crv.to_iges())
    _model.save(fn)
    if auto_view:
        view(fn)


'''1/4 Ellipse Arc'''
a = 10
b = 6
arc0 = ConicArc((0, -b, 0), (1, 0, 0), (a, 0, 0), (0, 1, 0), (a / math.sqrt(2), -b / math.sqrt(2), 0))
print(arc0(0), arc0(0, 1))
print(arc0(1), arc0(1, 1))
show(arc0, 'arc.igs')
