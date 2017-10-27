import numpy as np

from src.iges import IGES_Model
from src.nurbs.curve import ClampedNURBSCrv, GlobalInterpolatedCrv

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


pts = np.array([[0, 0, 0], [0, 3.14, 0], [0, 0, 4], [1, 0, -2],
                [0, 1, 0], [2, 0, 0], [0, 0, 9], [0.618, 1.414, 2.718]])

crv1 = GlobalInterpolatedCrv(pts, 3)
show(crv1, 'origin.igs')
print('Original curve:')
print(crv1.U)
print(crv1.weight)
print(crv1.cpt)

fn = 'split_test.igs'
model = IGES_Model(fn)

crv_list = ClampedNURBSCrv.split(crv1, [0.2, 0.6])
for k, crv in enumerate(crv_list):
    model.add_entity(crv.to_iges())
    print('Seg{}:'.format(k + 1))
    print(crv.U)
    print(crv.weight)
    print(crv.cpt)

model.write()
if auto_view:
    view(fn)
