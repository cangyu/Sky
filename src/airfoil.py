#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import math
import numpy as np
from numpy.linalg import norm
from matplotlib import pyplot as plt
from spacing import uniform, chebshev_dist
from misc import pnt_dist
from nurbs import GlobalInterpolatedCrv
from nurbs import RuledSurf
from settings import AIRFOIL_LIST, AIRFOIL_DIR

"""
Implementation of the Airfoil utilities.

Ref:
[1] https://github.com/dgorissen/naca
[2] https://github.com/adeharo9/NACA-airfoil-generator
[3] http://airfoiltools.com/airfoil/naca4digit
[4] http://airfoiltools.com/airfoil/naca5digit
"""


class AirfoilNotAvailable(Exception):
    pass


def naca4_parser(digits):
    m = float(digits[0]) / 100.0
    p = float(digits[1]) / 10.0
    t = float(digits[2:]) / 100.0
    return m, p, t


def naca4(m, p, t, n, trailing_blunt, half_cosine_spacing):
    """
    Calculate NACA 4-digit airfoil points.
    :param m: Max Camber.
    :type m: float
    :param p: Max Camber Position.
    :type p: float
    :param t: Relative Thickness.
    :type t: float
    :param n: Num of samples in X-direction.
    :type n: int
    :param trailing_blunt: Indicate if the trailing edge is closed or not.
    :type trailing_blunt: bool
    :param half_cosine_spacing: Indicate the distribution of sampling points.
    :type half_cosine_spacing: bool
    :return: Airfoil points(totally 2n+1) in XFOIL-style order.
    """

    def yc(x):
        if x < p:
            return m / p ** 2 * (2 * p * x - x ** 2)
        else:
            return m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x - x ** 2)

    def dyc(x):
        if x < p:
            return 2 * m / p ** 2 * (p - x)
        else:
            return 2 * m / (1 - p) ** 2 * (p - x)

    a0 = +0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = +0.2843
    a4 = -0.1015 if trailing_blunt else -0.1036

    def yt(x):
        return t / 0.2 * (a0 * x ** 0.5 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + a4 * x ** 4)

    def theta(x):
        return math.atan(dyc(x))

    def xu(x):
        return x - yt(x) * math.sin(theta(x))

    def yu(x):
        return yc(x) + yt(x) * math.cos(theta(x))

    def xl(x):
        return x + yt(x) * math.sin(theta(x))

    def yl(x):
        return yc(x) - yt(x) * math.cos(theta(x))

    u = chebshev_dist(0, 1, n + 1) if half_cosine_spacing else uniform(n + 1)

    upper_pnt = [(xu(_u), yu(_u)) for _u in u]
    lower_pnt = [(xl(_u), yl(_u)) for _u in u]

    return np.copy(upper_pnt[::-1] + lower_pnt[1:])


def naca5_parser(digits):
    cl = int(digits[0]) * (3.0 / 20.0)
    p = int(digits[1]) / 20.0
    q = bool(int(digits[2]))
    t = int(digits[3:]) / 100.0
    return cl, p, q, t


def naca5(cl, p, q, t, n, trailing_blunt, half_cosine_spacing):
    """
    Calculate NACA 5-digit airfoil points.
    :param cl: Designed coefficient of lift.
    :type cl: float
    :param p: Position of maximum camber.
    :type p: float
    :param q: Camber line reflex or not.
    :type q: bool
    :param t: Maximum thickness.
    :type t: float
    :param n: Num of samples in X-direction.
    :type n: int
    :param trailing_blunt: Indicate if the trailing edge is closed or not.
    :type trailing_blunt: bool
    :param half_cosine_spacing: Indicate the distribution of sampling points.
    :type half_cosine_spacing: bool
    :return: Airfoil points(totally 2n+1) in XFOIL-style order.
    """

    if not q:
        r = 3.33333333333212 * pow(p, 3) + 0.700000000000909 * pow(p, 2) + 1.19666666666638 * p - 0.00399999999996247
        k1 = 1514933.33335235 * pow(p, 4) - 1087744.00001147 * pow(p, 3) + 286455.266669048 * pow(p, 2) - 32968.4700001967 * p + 1420.18500000524

        def yc(x):
            ret = k1 / 6 * (pow(x, 3) - 3 * r * pow(x, 2) + pow(r, 2) * (3 - r) * x) if x < r else k1 / 6 * pow(r, 3) * (1 - x)
            return (cl / 0.3) * ret

        def dyc(x):
            ret = k1 / 6 * (3 * pow(x, 2) - 6 * r * x + pow(r, 2) * (3 - r)) if x < r else -k1 / 6 * pow(r, 3)
            return (cl / 0.3) * ret
    else:
        r = 10.6666666666861 * pow(p, 3) - 2.00000000001601 * pow(p, 2) + 1.73333333333684 * p - 0.0340000000002413
        k1 = -27973.3333333385 * pow(p, 3) + 17972.8000000027 * pow(p, 2) - 3888.40666666711 * p + 289.076000000022
        k21 = 85.5279999999984 * pow(p, 3) - 34.9828000000004 * pow(p, 2) + 4.80324000000028 * p - 0.21526000000003

        def yc(x):
            ret = k1 / 6 * (pow(x - r, 3) - k21 * pow(1 - r, 3) * x - pow(r, 3) * x + pow(r, 3)) if x < r else k1 / 6 * (k21 * pow(x - r, 3) - k21 * pow(1 - r, 3) * x + pow(r, 3) * (1 - x))
            return (cl / 0.3) * ret

        def dyc(x):
            ret = k1 / 6 * (3 * pow(x - r, 2) - k21 * pow(1 - r, 3) - pow(r, 3)) if x < r else k1 / 6 * (3 * k21 * pow(x - r, 2) - k21 * pow(1 - r, 3) - pow(r, 3))
            return (cl / 0.3) * ret

    a0 = +0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = +0.2843
    a4 = -0.1015 if trailing_blunt else -0.1036

    def yt(x):
        return t / 0.2 * (a0 * x ** 0.5 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + a4 * x ** 4)

    def theta(x):
        return math.atan(dyc(x))

    def xu(x):
        return x - yt(x) * math.sin(theta(x))

    def yu(x):
        return yc(x) + yt(x) * math.cos(theta(x))

    def xl(x):
        return x + yt(x) * math.sin(theta(x))

    def yl(x):
        return yc(x) - yt(x) * math.cos(theta(x))

    u = chebshev_dist(0, 1, n + 1) if half_cosine_spacing else uniform(n + 1)

    upper_pnt = [(xu(_u), yu(_u)) for _u in u]
    lower_pnt = [(xl(_u), yl(_u)) for _u in u]

    return np.copy(upper_pnt[::-1] + lower_pnt[1:])


class Airfoil(object):
    DEFAULT_AIRFOIL_NAME = 'Not Set'

    def __init__(self, pts, name=DEFAULT_AIRFOIL_NAME):
        """
        2D Airfoil, chord length is 1.
        :param pts:
        :param name:
        """

        self.name = name
        self.pts = np.copy(pts)

    @classmethod
    def from_file(cls, path, name=DEFAULT_AIRFOIL_NAME):
        p = np.loadtxt(path)
        if name == Airfoil.DEFAULT_AIRFOIL_NAME:
            name = os.path.basename(path)
            if name.endswith('.dat') or name.endswith('.txt'):
                base, ext = name.split('.')
                name = base
        return cls(p, name)

    @classmethod
    def from_local(cls, name):
        name = name.upper()
        if name in AIRFOIL_LIST:
            path = os.path.join(AIRFOIL_DIR, name + '.dat')
            p = np.loadtxt(path)
            return cls(p, name)
        else:
            raise AirfoilNotAvailable('{} not exist'.format(name))

    @classmethod
    def from_naca(cls, digits, n, trailing_blunt=True, half_cosine_spacing=True):
        """
        Calculate NACA 4/5-digit series airfoil points.
        :param digits: NACA airfoil specifier.
        :type digits: str
        :param n: Num of points.
        :type n: int
        :param trailing_blunt: Blunt flag.
        :param half_cosine_spacing: Spacing flag.
        :return: Coordinates assembled in the Airfoil object.
        :rtype: Airfoil
        """

        '''Check parameters and Conversions'''
        if digits.startswith('naca') or digits.startswith('NACA'):
            digits = digits[4:]
        assert digits.isdigit()

        n = n // 2 + 1

        trailing_blunt = bool(trailing_blunt)
        half_cosine_spacing = bool(half_cosine_spacing)

        '''Calculate points'''
        if len(digits) == 4:
            m, p, t = naca4_parser(digits)
            pts = naca4(m, p, t, n, trailing_blunt, half_cosine_spacing)
            return cls(pts, 'NACA' + digits)
        elif len(digits) == 5:
            l, p, q, t = naca5_parser(digits)
            pts = naca5(l, p, q, t, n, trailing_blunt, half_cosine_spacing)
            return cls(pts, 'NACA' + digits)
        else:
            raise AirfoilNotAvailable('NACA{} not resolvable.'.format(digits))

    def __repr__(self):
        return '{} with {} points'.format(self.name, self.size)

    @property
    def size(self):
        return len(self.pts)

    @property
    def crv(self):
        # return Spline(self.pts)
        return GlobalInterpolatedCrv(self.pts, 3)

    @property
    def tail_up(self):
        return self.pts[0]

    @property
    def tail_down(self):
        return self.pts[-1]

    @property
    def tail(self):
        return (self.tail_up + self.tail_down) / 2

    @property
    def front(self):
        """
        The most front point of the airfoil.
        """

        total = self.size
        cx = self.pts[0][0]
        k = 1
        while k < total and self.pts[k][0] < cx:
            cx = self.pts[k][0]
            k += 1

        return self.pts[k - 1]

    @property
    def chord_len(self):
        return pnt_dist(self.front, self.tail)

    @property
    def thickness(self):
        # TODO
        return 0.12

    @property
    def max_height(self):
        return self.thickness * self.chord_len

    @property
    def is_blunt(self):
        return not math.isclose(norm(self.tail_up - self.tail_down), 0)

    def save(self, fn):
        """
        Save all coordinates into file.
        :param fn: File name.
        :type fn: str
        :return: None.
        """

        f_out = open(fn, 'w')
        for p in self.pts:
            f_out.write('{:.8f}\t{:.8f}\n'.format(p[0], p[1]))
        f_out.close()

    def plot(self, ax):
        (px, py) = zip(*self.pts)
        ax.plot(px, py, '.-')
        ax.set_aspect('equal')

    def show(self):
        (px, py) = zip(*self.pts)
        plt.plot(px, py, '.-')
        plt.gca().set_aspect('equal')
        plt.show()

    def curvature_at(self, rel_pos):
        """
        Calculate the curvature at given position.
        :param rel_pos: Relative position.
        :type rel_pos: float
        :return: Curvature.
        :rtype: float
        """

        return self.curve.curvature(rel_pos)

    def refine(self, rel_pos):
        self.pts = self.crv.scatter(rel_pos)


def airfoil_interp(left_foil, right_foil, intermediate_pos, sample_pos):
    """
    Interpolate airfoils linearly.
    :param left_foil: Starting airfoil.
    :type left_foil: Airfoil
    :param right_foil: Ending airfoil.
    :type right_foil: Airfoil
    :param intermediate_pos: Relative positions between.
    :param sample_pos: Sampling positions on each curve.
    :return: Intermediate airfoils.
    """

    crv1 = left_foil.crv
    crv2 = right_foil.crv
    crv2.pan((0, 0, 1))
    rsf = RuledSurf(crv1, crv2)

    if type(intermediate_pos) in (float, int):
        crv = rsf.extract('v', intermediate_pos)
        assert np.greater_equal(sample_pos, 0) and np.less_equal(sample_pos, 1)
        pts = crv.scatter(sample_pos)
        return Airfoil(pts)
    elif type(intermediate_pos) in (np.ndarray, list):
        assert len(intermediate_pos) == len(sample_pos)
        ret = []
        for k, u in enumerate(intermediate_pos):
            crv = rsf.extract('v', u)
            csp = sample_pos[k]
            assert np.greater_equal(csp, 0).all() and np.less_equal(csp, 1).all()
            pts = crv.scatter(csp)
            ret.append(Airfoil(pts))
        return ret
    else:
        raise AssertionError('invalid input')


if __name__ == '__main__':
    naca0012 = Airfoil.from_naca('0012', 161)
    naca0012.save('NACA0012.dat')
    naca0012.show()

    naca23015 = Airfoil.from_naca('23015', 201)
    naca23015.save('NACA23015.dat')
    naca23015.show()

    naca23118 = Airfoil.from_naca('23118', 201)
    naca23118.save('NACA23118.dat')
    naca23118.show()

    naca13015 = Airfoil.from_naca('13015', 201)
    naca13015.save('NACA13015.dat')
    naca13015.show()
