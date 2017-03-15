import sys
import os
import math
import time
import numpy as np
import pylab as pl
from scipy import interpolate
from abc import ABCMeta, abstractmethod
from src.iges.iges_core import IGES_Model
from src.iges.iges_entity112 import *


z = np.array([0., 0.624029, 1.38967, 2.43503, 3.73439, 5.25574, 6.96162,
              8.81003, 10.7555, 12.75, 14.7445, 16.69, 18.5384, 20.2443, 21.7656,
              23.065, 24.1103, 24.876, 25.343, 25.5])

x_front = np.array([0., 0.05, 0.3, 1.7, 4.9, 6.85, 8.45, 9.65, 10.6, 11.1, 11.7, 12.1,
                    12.4, 12.8, 13.2, 13.7, 14.1, 14.5, 15.2, 16.])

x_tail = np.array([19.7, 19.6, 19.6, 19.5, 19.3, 19, 18.3, 17.3, 16.6, 16.5,
                   16.8, 17, 17.45, 17.8, 18.1, 18.4, 18.55, 18.65, 18.3, 17.8])

pl.plot(z, x_front, 'ro')
pl.plot(z, x_tail, 'bo')

z_new = np.linspace(0, 25.5, 1000)

xf_natural = interpolate.make_interp_spline(z, x_front, k=3, bc_type=([(2, 0)], [(2, 0)]))
xt_natural = interpolate.make_interp_spline(z, x_tail, k=3, bc_type=([(2, 0)], [(2, 0)]))

Coef = np.zeros((len(z), 3, 4))

for i in range(0, len(z)):
    Coef[i][2][0] = z[i]
    Coef[i][2][1] = float(1)

    Coef[i][0][0] = float(xf_natural(z[i]))
    Coef[i][0][1] = float(xf_natural(z[i], 1))
    Coef[i][0][2] = float(xf_natural(z[i], 2) / 2)
    Coef[i][0][3] = float(xf_natural(z[i], 3) / 6)

print(Coef)

front_curve = IGES_Entity112(z, Coef)
front_curve_record = front_curve.ConstructRecord()

xfn_new = xf_natural(z_new)
pl.plot(z_new, xfn_new, label="Front curve\n(Natural BC)")

xtn_new = xt_natural(z_new)
pl.plot(z_new, xtn_new, label="Tail curve\n(Natural BC)")

pl.legend(loc="lower right")
pl.show()
