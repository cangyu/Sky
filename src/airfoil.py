#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import math
import numpy as np
from numpy.linalg import norm
from matplotlib import pyplot as plt
import re
import collections
from spacing import uniform, chebshev_dist
from nurbs import Spline
from nurbs import RuledSurf
from settings import AIRFOIL_LIST, AIRFOIL_DIR, XFOIL_PATH

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


naca4_pattern = re.compile(r'^NACA(\d)(\d)(\d{2})$', re.IGNORECASE)


def naca4_parser(naca: str):
    param = re.search(naca4_pattern, naca).groups()
    m = int(param[0]) / 100.0
    p = int(param[1]) / 10.0
    t = int(param[2]) / 100.0
    return m, p, t


def naca4(m, p, t, n, trailing_blunt=True, half_cosine_spacing=True):
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
        return t / 0.2 * (a0 * pow(x, 0.5) + x * (a1 + x * (a2 + x * (a3 + x * a4))))

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

    u = chebshev_dist(0, 1, n) if half_cosine_spacing else uniform(n)

    upper_pnt = [(xu(_u), yu(_u)) for _u in u]
    lower_pnt = [(xl(_u), yl(_u)) for _u in u]

    return np.copy(upper_pnt[::-1] + lower_pnt[1:])


naca5_pattern = re.compile('^NACA(\d)(\d)(\d)(\d{2})$', re.IGNORECASE)


def naca5_parser(naca: str):
    param = re.search(naca5_pattern, naca).groups()
    cl = int(param[0]) * (3.0 / 20.0)
    p = int(param[1]) / 20.0
    q = bool(param[2])
    t = int(param[3]) / 100.0
    return cl, p, q, t


def naca5(cl, p, q, t, n, trailing_blunt=True, half_cosine_spacing=True):
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
        return t / 0.2 * (a0 * pow(x, 0.5) + x * (a1 + x * (a2 + x * (a3 + x * a4))))

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

    u = chebshev_dist(0, 1, n) if half_cosine_spacing else uniform(n)

    upper_pnt = [(xu(_u), yu(_u)) for _u in u]
    lower_pnt = [(xl(_u), yl(_u)) for _u in u]

    return np.copy(upper_pnt[::-1] + lower_pnt[1:])


naca6_pattern = re.compile(r'^NACA6(\d)\((\d)\)-(\d)(\d{2})$', re.IGNORECASE)


def naca6_parser(naca: str):
    param = re.search(naca6_pattern, naca).groups()
    series = int(param[0])
    low_drag_range = int(param[1]) / 10.0
    cl_design = int(param[2]) / 10.0
    thickness = int(param[3]) / 100.0
    return series, low_drag_range, cl_design, thickness


# TODO
def naca6(ser, rg, cl, t, n, a=1.0, half_cosine_spacing=True):
    g = -1.0 / (1 - a) * (pow(a, 2) * (0.5 * math.log(a) - 0.25) + 0.25)
    h = 1.0 / (1 - a) * (0.5 * pow(1 - a, 2) * math.log(1 - a) - 0.25 * (1 - pow(a, 2))) + g
    alpha = cl * h / (2 * math.pi * (a + 1))

    def yc(x):
        pass

    def dyc(x):
        pass

    def yt(x):
        pass


class Airfoil(object):
    DEFAULT_AIRFOIL_NAME = 'Not Set'

    def __init__(self, pts, name=DEFAULT_AIRFOIL_NAME):
        """
        2D Airfoil, chord length is 1.
        :param pts: Coordinates describing the airfoil.
        :param name: Name of the airfoil.
        :type name: str
        """

        self.name = name
        self.pts = np.copy(pts)

    @classmethod
    def from_file(cls, path, name=DEFAULT_AIRFOIL_NAME):
        if name == Airfoil.DEFAULT_AIRFOIL_NAME:
            name = os.path.basename(path)
            if name.endswith('.dat') or name.endswith('.txt'):
                base, ext = name.split('.')
                name = base

        p = np.loadtxt(path)
        return cls(p[:, :2], name)

    @classmethod
    def from_local(cls, name):
        name = name.upper()
        if name not in AIRFOIL_LIST:
            raise AirfoilNotAvailable('{} not exist'.format(name))

        path = os.path.join(AIRFOIL_DIR, name + '.dat')
        p = np.loadtxt(path)
        return cls(p[:, :2], name)

    @classmethod
    def from_naca(cls, digits, n, trailing_blunt=True, half_cosine_spacing=True):
        """
        Calculate NACA 4/5-digit series airfoil points.
        :param digits: NACA airfoil specifier.
        :type digits: str
        :param n: Num of points.
        :type n: int
        :param trailing_blunt: Blunt flag.
        :type trailing_blunt: bool
        :param half_cosine_spacing: Spacing flag.
        :type half_cosine_spacing: bool
        :return: Coordinates assembled in the Airfoil object.
        :rtype: Airfoil
        """

        n = n // 2 + 1
        foil = 'NACA' + digits

        '''Calculate points'''
        if re.match(naca4_pattern, foil) is not None:
            m, p, t = naca4_parser(foil)
            pts = naca4(m, p, t, n, trailing_blunt, half_cosine_spacing)
            return cls(pts, foil)
        elif re.match(naca5_pattern, foil) is not None:
            l, p, q, t = naca5_parser(foil)
            pts = naca5(l, p, q, t, n, trailing_blunt, half_cosine_spacing)
            return cls(pts, foil)
        else:
            raise AirfoilNotAvailable('{} not resolvable.'.format(foil))

    def __repr__(self):
        return '{} with {} points'.format(self.name, self.size)

    @property
    def size(self):
        return len(self.pts)

    def generate_nurbs_crv(self):
        return Spline(np.array([[p[0], p[1], 0] for p in self.pts]))

    @property
    def trailing_up(self):
        return self.pts[0]

    @property
    def trailing_down(self):
        return self.pts[-1]

    @property
    def trailing(self):
        return (self.trailing_up + self.trailing_down) / 2

    @property
    def is_blunt(self):
        return not math.isclose(norm(self.trailing_up - self.trailing_down), 0)

    def save(self, fn='', with_name=False):
        """
        Save all coordinates into file.
        :param fn: File name.
        :type fn: str
        :return: None.
        """

        if fn == '':
            fn = self.name + '.dat'

        f_out = open(fn, 'w')
        if with_name:
            f_out.write(self.name + '\n')

        for p in self.pts:
            f_out.write('{:10.6f}\t{:10.6f}\n'.format(p[0], p[1]))
        f_out.close()

    def plot(self, ax):
        (px, py) = zip(*self.pts)
        ax.plot(px, py, '.-')
        ax.set_aspect('equal')

    def show(self):
        (px, py) = zip(*self.pts)
        plt.plot(px, py, '.-')
        plt.gca().set_aspect('equal')
        plt.title(self.name)
        plt.show()

    def curvature_at(self, rel_pos):
        """
        Calculate the curvature at given position.
        :param rel_pos: Relative position.
        :return: Curvature.
        """

        crv = self.generate_nurbs_crv()

        if isinstance(rel_pos, collections.Iterable):
            return np.array([crv.curvature(u) for u in rel_pos])
        else:
            return crv.curvature(rel_pos)

    def refine(self, rel_pos):
        assert isinstance(rel_pos, collections.Iterable)
        crv = self.generate_nurbs_crv()
        self.pts = crv.scatter(rel_pos)[:, :2]


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

    crv1 = left_foil.generate_nurbs_crv()
    crv2 = right_foil.generate_nurbs_crv()
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


def find_alpha(foil, reynolds, ma, cl, iter_cnt=5000, alfa_only=True):
    """
    Find the AoA of given airfoil under specific conditions.
    :param foil: Target foil.
    :type foil: Airfoil
    :param reynolds: Reynolds Number.
    :type reynolds: float
    :param ma: Mach Number.
    :type ma: float
    :param cl: Target Lift Coefficient.
    :type cl: float
    :param iter_cnt: Num of iteration.
    :type iter_cnt: int
    :return: None.
    """

    foil_name = '_foil.dat'
    foil.save(foil_name)
    foil_path = os.path.join(os.getcwd(), foil_name)

    polar_fn = '_polar.dat'
    polar_path = os.path.join(os.getcwd(), polar_fn)

    cmd_fn = '_command.in'
    cmd_path = os.path.join(os.getcwd(), cmd_fn)
    cmd_stream = open(cmd_path, 'w')

    tmp_fn = '_tmp.out'
    tmp_path = os.path.join(os.getcwd(), tmp_fn)

    def issue_cmd(cmd):
        cmd_stream.write(cmd + '\n')

    '''Generate command file'''
    issue_cmd('load ' + foil_path)
    issue_cmd(foil.name)
    issue_cmd('panel')
    issue_cmd('oper')
    issue_cmd('visc {}'.format(reynolds))
    issue_cmd('M {}'.format(ma))
    issue_cmd('type 1')
    issue_cmd('pacc')
    issue_cmd(polar_fn)
    issue_cmd('')
    issue_cmd('iter')
    issue_cmd(str(iter_cnt))
    issue_cmd('cl {}'.format(cl))
    issue_cmd('')
    issue_cmd('')
    issue_cmd('quit')
    cmd_stream.close()

    '''Execute XFOIL commands'''
    os.system(XFOIL_PATH + ' < ' + cmd_path + ' > ' + tmp_path)
    os.remove(cmd_path)
    os.remove(foil_path)
    os.remove(tmp_path)

    '''Extract results'''
    polar_stream = open(polar_path)
    lines = polar_stream.readlines()
    polar_stream.close()
    os.remove(polar_path)

    '''Output'''
    results = lines[-1].split()
    alpha = float(results[0])
    foil_cl = float(results[1])
    foil_cd = float(results[2])
    foil_cm = float(results[4])
    foil_ld = foil_cl / foil_cd

    if alfa_only:
        return alpha
    else:
        return alpha, foil_cl, foil_cd, foil_cm, foil_ld


if __name__ == '__main__':
    naca0012 = Airfoil.from_naca('0012', 201)
    # naca0012.save()
    # naca0012.show()

    naca1408 = Airfoil.from_naca('1408', 131)
    # naca1408.save()
    # naca1408.show()

    naca13015 = Airfoil.from_naca('13015', 201)
    # naca13015.save()
    # naca13015.show()

    naca23015 = Airfoil.from_naca('23015', 201)
    # naca23015.save()
    # naca23015.show()

    naca23118 = Airfoil.from_naca('23118', 201)
    # naca23118.save()
    # naca23118.show()

    sc0712 = Airfoil.from_local('SC(2)-0712')
    # sc0712.show()
    aoa = find_alpha(sc0712, 6e7, 0.8, 0.91)
    print(aoa)
