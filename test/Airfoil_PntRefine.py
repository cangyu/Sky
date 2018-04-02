from matplotlib import pyplot as plt
from aircraft.wing import Airfoil
from grid import uniform

if __name__ == '__main__':
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
