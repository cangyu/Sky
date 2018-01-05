from matplotlib import pyplot as plt
from wing import Airfoil, airfoil_interp
from grid import uniform

if __name__ == '__main__':
    root_airfoil = Airfoil('NACA64(3)-218_160')
    tip_airfoil = Airfoil('NACA63(2)-615_161')
    mid = airfoil_interp(root_airfoil, tip_airfoil, [1 / 3, 2 / 3])
    n = 161
    nsp = uniform(n)
    p1 = mid[0].scatter(nsp)
    p2 = mid[1].scatter(nsp)

    foil1 = Airfoil(p1, name='mid1')
    foil1.save('interp1.dat')
    foil2 = Airfoil(p2)
    foil2.save('interp2.dat')

    fig = plt.figure()
    fig.set_size_inches(18.5, 16.5)
    root_ax = fig.add_subplot(411)
    root_airfoil.plot(root_ax)
    mid1_ax = fig.add_subplot(412, sharex=root_ax)
    foil1.plot(mid1_ax)
    mid2_ax = fig.add_subplot(413, sharex=root_ax)
    foil2.plot(mid2_ax)
    tip_ax = fig.add_subplot(414, sharex=root_ax)
    tip_airfoil.plot(tip_ax)
    fig.tight_layout()
    fig.savefig('test_airfoil_interp.png', dpi=300)
