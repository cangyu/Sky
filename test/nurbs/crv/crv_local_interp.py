from math import tan, sin, cos

import numpy as np

from src.iges import IGES_Model
from src.geom.curve import LocalCubicInterpolatedCrv

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

if __name__ == '__main__':
    cr = 12
    spn = 21

    '''Front'''
    fsl = np.array([1.2, 2.3, 2.45, 0])
    fsl[-1] = spn - sum(fsl[:-1])
    alpha = np.radians([50, 45, 33, 28])
    fp = np.zeros((5, 3))
    for i in range(4):
        fp[i + 1] = fp[i]
        fp[i + 1][0] += fsl[i] * tan(alpha[i])
        fp[i + 1][2] += fsl[i]
    ftv = np.array([[0, 0, 1],
                    [sin(alpha[1]), 0, cos(alpha[1])],
                    [sin(alpha[1]), 0, cos(alpha[1])],
                    [sin(alpha[3]), 0, cos(alpha[3])],
                    [sin(alpha[3]), 0, cos(alpha[3])]])

    '''Tail'''
    tsl = np.array([1.6, 4.0, 1.8, 0])
    tsl[-1] = spn - sum(tsl[:-1])
    beta = np.radians([-15, -48, -20, 15])
    tp = np.zeros((5, 3))
    tp[0][0] = cr
    for i in range(4):
        tp[i + 1] = tp[i]
        tp[i + 1][0] += tsl[i] * tan(beta[i])
        tp[i + 1][2] += tsl[i]
    ttv = np.array([[0, 0, 1],
                    [sin(beta[1]), 0, cos(beta[1])],
                    [sin(beta[1]), 0, cos(beta[1])],
                    [sin(beta[3]), 0, cos(beta[3])],
                    [sin(beta[3]), 0, cos(beta[3])]])

    '''Display'''
    fn = 'lci.igs'
    model = IGES_Model(fn)
    fc = LocalCubicInterpolatedCrv(fp, ftv)
    tc = LocalCubicInterpolatedCrv(tp, ttv)
    model.add_entity(fc.to_iges())
    model.add_entity(tc.to_iges())
    model.write()

    if auto_view:
        view(fn)
