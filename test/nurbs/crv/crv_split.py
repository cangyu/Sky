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

crv2 = crv1.decompose()

print(crv2.m)
print(crv2.U)
print(crv2.n)
print(crv2.Pw)
