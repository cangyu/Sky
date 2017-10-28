from copy import deepcopy

import numpy as np

from src.iges import IGES_Model
from src.geom.curve import ClampedNURBSCrv

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

    model = IGES_Model(fn)
    model.add_entity(crv.to_iges())
    model.write()
    if auto_view:
        view(fn)


U = np.array([0, 0, 0, 0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0])
Pw = np.array([[0, 3.14, 0, 1],
               [1, 0, -2, 1],
               [0, 1, 0, 1],
               [2, 0, 0, 1],
               [0, 0, 9, 1],
               [0.618, 1.414, 2.718, 1]])

C0 = ClampedNURBSCrv(U, Pw)
show(C0, 'C0.igs')
print('C0:')
print(C0.U)
print(C0.weight)

C1 = deepcopy(C0)
C1.reparameterization(2, 1, 3, 2)
show(C1, 'C1.igs')
print('C1:')
print(C0.U)
print(C0.weight)

C2 = deepcopy(C0)
C2.reparameterization(2, 0, 3, 1)
show(C2, 'C2.igs')
print('C2:')
print(C2.U)
print(C2.weight)
