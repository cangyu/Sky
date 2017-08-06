from src.nurbs.curve import ConicArc
from src.iges.iges_core import IGES_Model

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
    :type crv: ClampedNURBSCrv
    :param fn: File name
    :type fn: str
    :return: None
    """

    _model = IGES_Model(fn)
    _model.add_entity(crv.to_iges())
    _model.write()
    if auto_view:
        view(fn)


arc0 = ConicArc((0, 0, 0), (0, 1, 0), (4.66, 2.41, 0), (1, 0, 0), (2.2, 2.2, 0))
show(arc0, 'arc0.igs')
