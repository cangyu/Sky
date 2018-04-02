from matplotlib import pyplot as plt
import numpy as np
from aircraft.wing import WingProfile, Airfoil

if __name__ == '__main__':
    af = Airfoil('NACA0012')
    ending = np.array([(15, 10, 0), (0, 0, 0)])
    wp = WingProfile(af, ending)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    wp.plot(ax)
    plt.show()
