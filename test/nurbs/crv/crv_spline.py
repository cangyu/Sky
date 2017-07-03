import unittest
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from src.nurbs.curve import Spline, GlobalInterpolatedCrv
from src.iges.iges_core import IGES_Model
from settings import AIRFOIL_DIR

try:
    from src.com.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = False


def get_airfoil_pts(airfoil):
    fn = AIRFOIL_DIR + '/' + airfoil + '.dat'
    fin = open(fn)
    pnt_list = []
    for pnt in fin:
        x, y, z = pnt.split()
        pnt_list.append([float(x), float(y), float(z)])
    fin.close()

    return np.copy(pnt_list)


def write_airfoil(airfoil):
    spline_foil = Spline(get_airfoil_pts(airfoil))
    fn = airfoil + '_spline' + '.igs'
    model_file = IGES_Model(fn)
    model_file.add_entity(spline_foil.to_iges())
    model_file.write()
    if auto_view:
        view(fn)


def compare_airfoil(airfoil, d, N=1000):
    pts = get_airfoil_pts(airfoil)
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
