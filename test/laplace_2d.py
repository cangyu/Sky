from src.msh.elliptical import Laplace_2D


def airfoil(foil, L, R):
    pass


class Laplace2D_Test(object):
    def test_airfoil(self):
        msh = airfoil("NACA0012", 10, 50)
        U, V = 60, 15
        show_uniform(msh, U, V)
        write_uniform_p3d(msh, U, V, "NACA0012.xyz")
