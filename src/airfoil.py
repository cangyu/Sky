#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import math
import numpy as np
from numpy.linalg import norm
from grid import LinearTFI2D, Laplace2D, ThomasMiddlecoff2D
from grid import Plot3D, Plot3DBlock
from spacing import uniform, chebshev_dist
from spacing import hyperbolic_tangent, single_exponential, double_exponential
from iges import Model, Entity116
from misc import pnt_dist
from nurbs import Line, ConicArc
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

    upper_pnt = np.array([(xu(_u), yu(_u)) for _u in u])
    lower_pnt = np.array([(xl(_u), yl(_u)) for _u in u])

    return upper_pnt[::-1] + lower_pnt[1:]


def naca5_parser(digits):
    cld = int(digits[0]) * (3.0 / 2.0) / 10.0
    p = 0.5 * int(digits[1]) / 100.0
    q = bool(int(digits[2]))
    t = int(digits[3:]) / 100.0
    return cld, p, q, t


def naca5(cl, p, q, t, n, trailing_blunt, half_cosine_spacing):
    """
    Calculate NACA 5-digit airfoil points.
    :param digits: Airfoil specifier.
    :type digits: str
    :param n: Num of samples in X-direction.
    :type n: int
    :param trailing_blunt: Indicate if the trailing edge is closed or not.
    :type trailing_blunt: bool
    :param half_cosine_spacing: Indicate the distribution of sampling points.
    :type half_cosine_spacing: bool
    :return: Airfoil points(totally 2n+1) in XFOIL-style order.
    """

    a0 = +0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = +0.2843
    a4 = -0.1015 if trailing_blunt else -0.1036

    if half_cosine_spacing:
        beta = linspace(0.0, pi, n + 1)
        x = [(0.5 * (1.0 - cos(x))) for x in beta]  # Half cosine based spacing
    else:
        x = linspace(0.0, 1.0, n + 1)

    yt = [5 * t * (a0 * sqrt(xx) + a1 * xx + a2 * pow(xx, 2) + a3 * pow(xx, 3) + a4 * pow(xx, 4)) for xx in x]

    P = [0.05, 0.1, 0.15, 0.2, 0.25]
    M = [0.0580, 0.1260, 0.2025, 0.2900, 0.3910]
    K = [361.4, 51.64, 15.957, 6.643, 3.230]

    m = interpolate(P, M, [p])[0]
    k1 = interpolate(M, K, [m])[0]

    xc1 = [xx for xx in x if xx <= p]
    xc2 = [xx for xx in x if xx > p]
    xc = xc1 + xc2

    if p == 0:
        xu = x
        yu = yt

        xl = x
        yl = [-x for x in yt]

        zc = [0] * len(xc)
    else:
        yc1 = [k1 / 6.0 * (pow(xx, 3) - 3 * m * pow(xx, 2) + pow(m, 2) * (3 - m) * xx) for xx in xc1]
        yc2 = [k1 / 6.0 * pow(m, 3) * (1 - xx) for xx in xc2]
        zc = [cld / 0.3 * xx for xx in yc1 + yc2]

        dyc1_dx = [cld / 0.3 * (1.0 / 6.0) * k1 * (3 * pow(xx, 2) - 6 * m * xx + pow(m, 2) * (3 - m)) for xx in xc1]
        dyc2_dx = [cld / 0.3 * (1.0 / 6.0) * k1 * pow(m, 3)] * len(xc2)

        dyc_dx = dyc1_dx + dyc2_dx
        theta = [atan(xx) for xx in dyc_dx]

        xu = [xx - yy * sin(zz) for xx, yy, zz in zip(x, yt, theta)]
        yu = [xx + yy * cos(zz) for xx, yy, zz in zip(zc, yt, theta)]

        xl = [xx + yy * sin(zz) for xx, yy, zz in zip(x, yt, theta)]
        yl = [xx - yy * cos(zz) for xx, yy, zz in zip(zc, yt, theta)]

    X = xu[::-1] + xl[1:]
    Z = yu[::-1] + yl[1:]

    return X, Z


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
        :param n: Num of points.
        :param trailing_blunt: Blunt flag.
        :param half_cosine_spacing: Spacing flag.
        :return: Coordinates assembled in the Airfoil object.
        :rtype: Airfoil
        """

        '''Check parameters and Conversions'''
        assert type(digits) is str
        if digits.startswith('naca') or digits.startswith('NACA'):
            digits = digits[4:]
        assert digits.isdigit()

        assert type(n) is int
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

        total = self.pnt_num
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

    def to_blunt(self):
        pfx = self.front[0]
        r0 = self.chord_len
        self.pts = self.pts[1:-1]
        r1 = self.chord_len
        ratio = r0 / r1
        for k in range(len(self.pts)):
            self.pts[k][0] = pfx + (self.pts[k][0] - pfx) * ratio

    def save(self, fn):
        """
        Save all coordinates into file.
        :param fn: File name.
        :type fn: str
        :return: None.
        """

        f_out = open(fn, 'w')
        for p in self.pts:
            f_out.write('{:10.6f}\t{:10.6f}\t{:10.6f}\n'.format(p[0], p[1], p[2]))
        f_out.close()

    def plot(self, ax):
        (px, py, pz) = zip(*self.pts)
        ax.plot(px, py, '.-')
        ax.set_aspect('equal')

    def curvature_at(self, rel_pos):
        """
        Calculate the curvature at given position.
        :param rel_pos: Relative position.
        :type rel_pos: float
        :return: Curvature.
        :rtype: float
        """

        return self.curve.curvature(rel_pos)

    def gen_grid(self, *args, **kwargs):
        """
        Generate grid for 2D airfoil or wing profile.
        :param args: Containing the geometric description and node distribution of the grid.
        :param kwargs: Extra options on smoothing and spacing.
        :return: The wire-frame, plot3d-grid and fluent-grid(with predefined BC) of the flow field.
        """

        '''
        a: Width of the front part of the flow field.
        b: Semi-height of the flow field.
        c: Width of the rear part of the flow field.
        n0: Num of points along airfoil.
        n1: Num of points along vertical direction.
        n2: Num of points along horizontal direction in rear field.
        n3: Num of Points on the trailing edge.
        '''
        assert len(args) == 7
        a, b, c, n0, n1, n2, n3 = args

        wire_frame = Model()
        p3d_grid = Plot3D()

        '''Flow-field geometries'''
        pts = np.empty((8, 3), float)
        pts[0] = self.tail_up
        pts[1] = self.tail_down
        pts[2] = np.array([0, b, self.z])
        pts[3] = np.array([0, -b, self.z])
        pts[4] = np.array([c, pts[0][1], self.z])
        pts[5] = np.array([c, pts[1][1], self.z])
        pts[6] = np.array([c, b, self.z])
        pts[7] = np.array([c, -b, self.z])

        crv = [self.crv,  # c0
               ConicArc(pts[2], (-1, 0, 0), pts[3], (1, 0, 0), (-a, 0, 0)),  # c1
               Line(pts[0], pts[2]),  # c2
               Line(pts[1], pts[3]),  # c3
               Line(pts[4], pts[6]),  # c4
               Line(pts[5], pts[7]),  # c5
               Line(pts[0], pts[4]),  # c6
               Line(pts[1], pts[5]),  # c7
               Line(pts[2], pts[6]),  # c8
               Line(pts[3], pts[7]),  # c9
               Line(pts[0], pts[1]),  # c10
               Line(pts[4], pts[5])]  # c11

        '''Construct wire-frame'''
        for p in pts:
            wire_frame.add(Entity116(p[0], p[1], p[2]))
        for c in crv:
            wire_frame.add(c.to_iges())

        '''Knot distribution'''
        u = [double_exponential(n0, 0.5, -1.5, 0.5),  # c0, c1
             hyperbolic_tangent(n1, 2),  # c2, c3, c4, c5
             single_exponential(n2, 3),  # c6, c7, c8, c9
             uniform(n3)]  # c10, c11

        '''Structured grid blocks'''
        leading_blk = LinearTFI2D(crv[2], crv[0], crv[3], crv[1])
        tailing_up_blk = LinearTFI2D(crv[6], crv[2], crv[8], crv[4])
        tailing_down_blk = LinearTFI2D(crv[3], crv[7], crv[5], crv[9])
        rear_blk = LinearTFI2D(crv[10], crv[6], crv[11], crv[7])

        '''Construct Plot3D grid for basic checking'''
        leading_blk.calc_grid(u[1], u[0])
        leading_tfi_grid = leading_blk.grid
        leading_smooth_ok = False
        if 'leading_smooth' in kwargs:
            smooth = kwargs['leading_smooth']
            if smooth in ('Laplace', 'laplace'):
                leading_grid_laplace = Laplace2D(leading_tfi_grid)
                leading_grid_laplace.smooth()
                p3d_grid.add(Plot3DBlock.construct_from_array(leading_grid_laplace.grid))
                leading_smooth_ok = True
            if smooth in ('TM', 'tm', 'Thomas-Middlecoff', 'thomas-middlecoff'):
                leading_grid_tm = ThomasMiddlecoff2D(leading_tfi_grid)
                leading_grid_tm.smooth()
                p3d_grid.add(Plot3DBlock.construct_from_array(leading_grid_tm.grid))
                leading_smooth_ok = True
        if not leading_smooth_ok:
            p3d_grid.add(Plot3DBlock.construct_from_array(leading_tfi_grid))

        tailing_up_blk.calc_grid(u[2], u[1])
        tailing_up_grid = tailing_up_blk.grid
        p3d_grid.add(Plot3DBlock.construct_from_array(tailing_up_grid))

        tailing_down_blk.calc_grid(u[1], u[2])
        tailing_down_grid = tailing_down_blk.grid
        p3d_grid.add(Plot3DBlock.construct_from_array(tailing_down_grid))

        if self.is_blunt:
            rear_blk.calc_grid(u[3], u[2])
            rear_grid = rear_blk.grid
            p3d_grid.add(Plot3DBlock.construct_from_array(rear_grid))

        return wire_frame, p3d_grid

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
