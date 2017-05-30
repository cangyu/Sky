import unittest
import numpy as np
import math
from src.aircraft.general import GeneralFrame
from src.aircraft.wing import Wing


def chebshev_dist(start, end, n):
    ang = np.linspace(math.pi, 0, n)
    pr = np.zeros(n)
    for i in range(0, n):
        pr[i] = math.cos(ang[i])
        pr[i] = start + (end - start) / 2 * (pr[i] + 1)

    return pr


class bwb_section_surf_test(unittest.TestCase):
    def test(self):
        gf = GeneralFrame()

        airfoil = ['M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6', 'M6']
        thkf = np.ones(len(airfoil))
        z = chebshev_dist(0, gf.Bt, len(airfoil))
        xf = gf.xfront(z)
        yf = gf.yfront(z)
        xt = gf.xtail(z)
        yt = gf.ytail(z)

        wg = Wing(airfoil, thkf, z, xf, yf, xt, yt)
        wg.write('BWB.igs')
