import unittest
from matplotlib import pyplot as plt
from wing import Airfoil
from grid import uniform


class AirfoilTestCase(unittest.TestCase):
    def test_grid(self):
        # airfoil, A, B, C, N0, N1, N2, N3
        data = [('SC(2)-0406', 30, 20, 50, 90, 60, 80, 3, 'none'),
                ('RAE2822', 30, 20, 50, 90, 60, 80, 3, 'none'),
                ('SC(2)-0406', 30, 20, 50, 90, 60, 80, 3, 'laplace'),
                ('RAE2822', 30, 20, 50, 90, 60, 80, 3, 'laplace'),
                ('NLF(1)-0414F', 30, 20, 50, 91, 61, 80, 3, 'thomas-middlecoff'),
                ('RAE2822', 30, 20, 50, 90, 60, 80, 3, 'thomas-middlecoff')]

        for k in range(len(data)):
            fn, la, lb, lc, n0, n1, n2, n3, smt = data[k]
            foil = Airfoil(fn)
            bunch = foil.gen_grid(la, lb, lc, n0, n1, n2, n3, leading_smooth=smt)
            p3d = bunch[1]
            p3d.save(fn + '_flowfield_grid-smooth={}.xyz'.format(smt))
        self.assertTrue(True)

    def test_refine(self):
        airfoil = 'NACA64(3)-218'
        refined_num = 161
        foil_compare_fig = plt.figure()
        original_ax = foil_compare_fig.add_subplot(211)
        test_foil = Airfoil(airfoil)
        test_foil.plot(original_ax)
        original_ax.set_title(airfoil + ' original')
        new_sp = uniform(refined_num)
        test_foil.refine(new_sp)
        current_ax = foil_compare_fig.add_subplot(212)
        test_foil.plot(current_ax)
        current_ax.set_title(airfoil + ' refined')
        test_foil.save('{}_{}.dat'.format(airfoil, refined_num))
        foil_compare_fig.savefig('test_airfoil_refine.png', dpi=600)
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
