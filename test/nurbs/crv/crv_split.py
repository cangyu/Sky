import numpy as np
from src.nurbs.curve import ClampedNURBSCrv, GlobalInterpolatedCrv, Spline, calc_pnt_param
from src.iges.iges_core import IGES_Model

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

fn = 'test.igs'
model = IGES_Model(fn)

pts = np.array([[0, 0, 0],
                [0, 3.14, 0],
                [0, 0, 4],
                [1, 0, -2],
                [0, 1, 0],
                [2, 0, 0],
                [0, 0, 9],
                [0.618, 1.414, 2.718]])

crv1 = GlobalInterpolatedCrv(pts, 3)

crv_list = ClampedNURBSCrv.split(crv1, [0.2, 0.6])
for k, crv in enumerate(crv_list):
    print('Seg{}:'.format(k + 1))
    print(crv.U)
    print(crv.weight)
    print(crv.cpt)
