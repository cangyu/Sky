import unittest

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from iges import Model
from nurbs import Spline, GlobalInterpolatedCrv
from src.aircraft.wing import Airfoil

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True


def get_airfoil_pts(foil):
    return Airfoil.read_pts(foil)


def write_airfoil(airfoil):
    spline_foil = Spline(get_airfoil_pts(airfoil))
    fn = airfoil + '_spline' + '.igs'
    model_file = Model()
    model_file.add(spline_foil.to_iges())
    model_file.save(fn)
    if auto_view:
        view(fn)


def compare_airfoil(foil, d, N=1000):
    pts = get_airfoil_pts(foil)
    spline_foil = Spline(pts)
    ginterp_foil = GlobalInterpolatedCrv(pts)

    u = np.linspace(0, 1, N)
    v = np.empty((2, N, 3), float)
    for i in range(N):
        v[0][i] = spline_foil(u[i], d)
        v[1][i] = ginterp_foil(u[i], d)

    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(v[0, :, 0], v[0, :, 1], v[0, :, 2], label='Spline_{}'.format(d))
    ax.plot(v[1, :, 0], v[1, :, 1], v[1, :, 2], label='GlobalInterp_{}'.format(d))
    ax.legend()
    plt.show()


class SplineInterpTest(unittest.TestCase):
    airfoil_list = ['M6', 'NACA0012', 'RAE2822']
    pder_list = [0, 1, 2]

    @staticmethod
    def test_basic():
        print('Basic Testing:')
        for foil in SplineInterpTest.airfoil_list:
            print(foil)
            write_airfoil(foil)

    @staticmethod
    def test_compare():
        print('Compare with GlobalInterpolatedCrv:')
        for pder in SplineInterpTest.pder_list:
            print('pder:{}'.format(pder))
            for foil in SplineInterpTest.airfoil_list:
                print(foil)
                compare_airfoil(foil, pder)


if __name__ == '__main__':
    unittest.main()
