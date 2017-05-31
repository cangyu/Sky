import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import interpolate


class GeneralFrame(object):
    def __init__(self, param):
        self.Cr = float(param[0])
        self.Cm = float(param[1])
        self.Ct = float(param[2])
        self.Bm = float(param[3])
        self.Bt = float(param[4])
        self.Alpha_r = float(param[5])
        self.Alpha_m = float(param[6])
        self.Alpha_t = float(param[7])

        self.front_pnt = np.zeros((3, 3), float)
        self.tail_pnt = np.zeros((3, 3), float)

        self.xfront = None
        self.yfront = None
        self.xtail = None
        self.ytail = None

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
        yf = np.array([self.front_pnt[0][1], self.front_pnt[1][1], self.front_pnt[2][1]])
        xt = np.array([self.tail_pnt[0][0], self.tail_pnt[1][0], self.tail_pnt[2][0]])
        yt = np.array([self.tail_pnt[0][1], self.tail_pnt[1][1], self.tail_pnt[2][1]])
        self.xfront = interpolate.make_interp_spline(z, xf, 3, bc_type=([(1, 0)], [(2, 0)]))
        self.yfront = interpolate.make_interp_spline(z, yf, 3, bc_type=([(1, 0)], [(2, 0)]))
        self.xtail = interpolate.make_interp_spline(z, xt, 3, bc_type=([(1, 0)], [(2, 0)]))
        self.ytail = interpolate.make_interp_spline(z, yt, 3, bc_type=([(1, 0)], [(2, 0)]))

    def plot(self, N=1000):
        z = np.linspace(0, self.Bt, N)
        xf = self.xfront(z)
        xt = self.xtail(z)

        plt.figure()
        plt.plot(z, xf, label='Front')
        plt.plot(z, xt, label='Tail')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.show()
