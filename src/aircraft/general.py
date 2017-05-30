import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import interpolate


class GeneralFrame(object):
    def __init__(self):
        self.Cr = 100.0
        self.Cm = 60.0
        self.Ct = 20.0
        self.Bm = 30.0
        self.Bt = 105.0
        self.Alpha_r = 0
        self.Alpha_m = 45
        self.Alpha_t = 30

        self.front_pnt = np.zeros((3, 3))
        self.tail_pnt = np.zeros((3, 3))

        self.xfront = None
        self.xtail = None

        self.update()

    def update(self):
        self.front_pnt[0] = np.zeros(3)
        self.tail_pnt[0] = np.array([self.Cr, 0, 0])
        self.front_pnt[1] = np.array([self.Bm * math.tan(math.radians(self.Alpha_m)), 0, self.Bm])
        self.tail_pnt[1] = np.array([self.front_pnt[1][0] + self.Cm, 0, self.Bm])
        self.front_pnt[2] = np.array([self.front_pnt[1][0] + (self.Bt - self.Bm) * math.tan(math.radians(self.Alpha_t)), 0, self.Bt])
        self.tail_pnt[2] = np.array([self.front_pnt[2][0] + self.Ct, 0, self.Bt])

        z = np.array([0, self.Bm, self.Bt])
        xf = np.array([self.front_pnt[0][0], self.front_pnt[1][0], self.front_pnt[2][0]])
        xt = np.array([self.tail_pnt[0][0], self.tail_pnt[1][0], self.tail_pnt[2][0]])
        self.xfront = interpolate.make_interp_spline(z, xf, 3, bc_type=([(1, 0)], [(2, 0)]))
        self.xtail = interpolate.make_interp_spline(z, xt, 3, bc_type=([(1, 0)], [(2, 0)]))

    def plot(self):
        z = np.linspace(0, self.Bt, 1000)
        xf = self.xfront(z)
        xt = self.xtail(z)

        plt.figure()
        plt.plot(z, xf, label='Front')
        plt.plot(z, xt, label='Tail')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.show()
