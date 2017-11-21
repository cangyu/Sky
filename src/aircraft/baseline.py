import unittest
import math
from copy import deepcopy
import numpy as np
from abc import ABCMeta, abstractmethod
from math import sin, cos, tan, radians, fabs, atan2
from scipy.integrate import romberg
from scipy.optimize import root
from nurbs import Crv, ConicArc, to_homogeneous, to_cartesian, LocalCubicInterpolatedCrv, Line, point_inverse, Circle, Skinned, RuledSurf, Coons
from misc import sqrt2


class Nose(object):
    def __init__(self, nl, nfr, dh, h2, w2):
        """
        Simplified nose of a 'Wing-Tube' configuration aircraft.
        :param nl: Length of nose.
        :type nl: float
        :param nfr: Radius of the front hole.
        :type nfr: float
        :param dh:  Delta height between the center of the front profile and the back profile.
        :type dh: float
        :param h2: Half of the height of the back profile.
        :type h2: float
        :param w2: Half of the width of the back profile.
        :type w2: float
        """

        a = h2
        b = w2

        p0 = np.array([0, nfr, 0])
        t0 = np.array([0, 1, 0])
        p2 = np.array([nl, dh + a, 0])
        t2 = np.array([1, 0, 0])
        p1c = np.array([nl, nfr, 0])
        p1 = p1c + np.array([nl * math.cos(math.radians(135)), (p2[1] - nfr) * math.sin(math.radians(135)), 0])
        arc1 = ConicArc(p0, t0, p2, t2, p1)

        p3 = np.array([0, -nfr, 0])
        t3 = np.array([0, -1, 0])
        p5 = p2 - np.array([0, 2 * a, 0])
        t5 = np.array([1, 0, 0])
        p4c = np.array([nl, -nfr, 0])
        p4 = p4c + np.array([nl * math.cos(math.radians(225)), (-nfr - p5[1]) * math.sin(math.radians(225)), 0])
        arc2 = ConicArc(p3, t3, p5, t5, p4)

        p9c = (p2 + p5) / 2
        p6 = np.array([0, 0, nfr])
        t6 = (0, 0, 1)
        p9 = np.array([0, 0, b]) + p9c
        p8 = np.array([p9[0], 0, p9[2]])
        t8 = (1, 0, 0)
        p7 = np.array([nl, 0, nfr]) + np.array([-nl / sqrt2, 0, (b - nfr) / sqrt2])

        arc3 = Circle.from_2pnt(p0, p3, 180, (1, 0, 0))
        arc5 = ConicArc(p6, t6, p8, t8, p7)
        for i in range(3):
            arc5.Pw[i][1] = i / 2 * p9[1] * arc5.Pw[i][-1]
        arc5.reset(arc5.U, arc5.Pw)

        '''Construct the ellipse by stretching the circle'''
        arc4 = Circle.from_2pnt(p2, p5, 180, (1, 0, 0))
        arc4.reset(arc4.U, np.copy(list(map(lambda u: to_homogeneous(to_cartesian(u) * (1, 1, b / a), u[-1]), arc4.Pw))))

        arc3_1, arc3_2 = Crv.split(arc3, [0.5])
        arc4_1, arc4_2 = Crv.split(arc4, [0.5])

        self.frame_up = deepcopy(arc1)
        self.frame_down = deepcopy(arc2)
        self.frame_mid = deepcopy(arc5)
        self.frame_front_up = deepcopy(arc3_1)
        self.frame_front_down = deepcopy(arc3_2)
        self.frame_back_up = deepcopy(arc4_1)
        self.frame_back_down = deepcopy(arc4_2)


class Fuselage(object):
    def __init__(self, l1, l2, fl):
        """
        Simplified fuselage of a 'Wing-Tube' configuration aircraft.
        :param l1: Upper part of the front profile curve.
        :type l1: Crv
        :param l2: Lower part of the front profile curve.
        :type l2: Crv
        :param fl: Length of the fuselage.
        :type fl: float
        """

        arc6_1 = deepcopy(l1)
        arc6_2 = deepcopy(l2)
        arc6_1.pan((fl, 0, 0))
        arc6_2.pan((fl, 0, 0))

        line1 = Line(l1(0), arc6_1(0))
        line2 = Line(l2(0), arc6_2(0))
        line3 = Line(l2(1), arc6_2(1))

        self.frame_front_up = deepcopy(l1)
        self.frame_front_down = deepcopy(l2)
        self.frame_back_up = deepcopy(arc6_1)
        self.frame_back_down = deepcopy(arc6_2)
        self.frame_up = deepcopy(line1)
        self.frame_mid = deepcopy(line2)
        self.frame_down = deepcopy(line3)

        fuselage_surf_up = Coons(l1, arc6_1, line1, line2)
        self.surf_up = deepcopy(fuselage_surf_up)

        fuselage_surf_down = Coons(l2, arc6_2, line2, line3)
        self.surf_down = deepcopy(fuselage_surf_down)


class Tail(object):
    def __init__(self, fu, fd, tl, dh, ttr):
        """
        Simplified tail of a 'Wing-Tube' configuration aircraft.
        :param fu: Upper part of the front profile.
        :type fu: Crv
        :param fd: Lower part of the front profile.
        :type fd: Crv
        :param tl: Length of the tail.
        :type tl: float
        :param dh: Delta of the upper frame on back profile.
        :type dh: float
        :param ttr: Radius of the hole on back profile.
        :type ttr: float
        """

        self.frame_front_up = deepcopy(fu)
        self.frame_front_down = deepcopy(fd)

        p10 = fu.start
        t10 = (1, 0, 0)
        p11c = p10 - np.array([0, dh, 0])
        p12 = p11c + np.array([tl, 0, 0])
        t12 = (0, -1, 0)
        p11 = p11c + np.array([tl / sqrt2, dh / sqrt2, 0])
        arc6 = ConicArc(p10, t10, p12, t12, p11)
        self.frame_up = deepcopy(arc6)

        p13 = fd.end
        t13 = (1, 0, 0)
        p15 = p12 - np.array([0, 2 * ttr, 0])
        t15 = (0, 1, 0)
        p14c = np.array([p13[0], p15[1], 0])
        p14 = p14c + np.array([tl / sqrt2, -(p15[1] - p13[1]) / sqrt2, 0])
        arc7 = ConicArc(p13, t13, p15, t15, p14)
        self.frame_down = deepcopy(arc7)

        arc8 = Circle.from_2pnt(p12, p15, 180, (1, 0, 0))
        arc8_1, arc8_2 = ClampedNURBSCrv.split(arc8, [0.5])
        self.frame_back_up = deepcopy(arc8_1)
        self.frame_back_down = deepcopy(arc8_2)

        p16 = fd.start
        t16 = (1, 0, 0)
        p18 = arc8_2.start
        p19 = np.array([p18[0], p16[1], p18[2]])
        t19 = (0, 0, -1)
        p17c = np.array([p16[0], p16[1], p19[2]])
        p17 = p17c + np.array([tl / sqrt2, 0, (p16[2] - p19[2]) / sqrt2])

        arc9 = ConicArc(p16, t16, p19, t19, p17)
        tmp = p18[1] - p16[1]
        n = len(arc9.Pw)
        for k, pw in enumerate(arc9.Pw):
            w = pw[-1]
            a, b, c = to_cartesian(pw)
            b += k / (n - 1) * tmp
            arc9.Pw[k] = to_homogeneous((a, b, c), w)
        arc9.reset(arc9.U, arc9.Pw)

        self.frame_mid = deepcopy(arc9)


class WingPlanform(metaclass=ABCMeta):
    @abstractmethod
    def x_front(self, u):
        """
        Calculate the X-coordinate on leading edge.
        :param u: Relative position parameter.
        :type u: float
        :return: X-coordinate on leading edge.
        :rtype: float
        """

        pass

    @abstractmethod
    def x_tail(self, u):
        """
        Calculate the X-coordinate on trailing edge.
        :param u: Relative position parameter.
        :type u: float
        :return: X-coordinate on trailing edge.
        :rtype: float
        """

        pass

    @abstractmethod
    def z(self, u):
        """
        Calculate the Z-coordinate in span-wise direction.
        :param u: Relative position parameter.
        :type u: float
        :return: Z-coordinate in span-wise direction.
        :rtype: float
        """

        pass

    def chord_len(self, u):
        return self.x_tail(u) - self.x_front(u)

    @property
    def area(self):
        """
        Total area of the planar wing. (Not a half)
        :return: The area.
        :rtype: float
        """

        return 2 * fabs(romberg(self.chord_len, 0, 1)) * self.z(1)

    @property
    def mean_aerodynamic_chord(self):
        """
        Get the mean aerodynamic chord length of the wing.
        :return: The MAC.
        :rtype: float
        """

        return self.z(1) * fabs(romberg(lambda u: self.chord_len(u) ** 2, 0, 1)) / (0.5 * self.area)

    @property
    def span(self):
        return 2 * (self.z(1) - self.z(0))

    @property
    def root_chord_len(self):
        return self.chord_len(0)

    @property
    def tip_chord_len(self):
        return self.chord_len(1)


    @abstractmethod
    def __repr__(self):
        pass


class HorizontalStablizer(object):
    def __init__(self):
        self.surf = None
        self.tail = None

    @classmethod
    def from_section_param(cls, foil, zoff, cl, swp_bk, tws, dih, ptws, rfy, tkf, pan_dir):
        """
        Construct the vertical stablizer from profile geometrical parameters.
        :param foil: Section airfoil list.
        :param zoff: Section 'z'-dim offset list when constructing from wing.
        :param cl: Section chord length list.
        :param swp_bk: Section sweep back angle list.
        :param tws:  Section twist angle list.
        :param dih: Section dihedral angle list.
        :param ptws: Section twist relative position list.
        :param rfy: Section 'y'-dim ref list when constructing from wing.
        :param tkf: Section thickness factor list.
        :param pan_dir: Panning direction.
        :return: The vertical stablizer.
        :rtype: HorizontalStablizer
        """

        '''Defensive check'''
        n = len(foil)
        if n < 2:
            raise ValueError("Insufficient input.")

        crv_list = []
        for k in range(n):
            wp = WingProfile.from_geom_param(foil[k], zoff[k], cl[k], swp_bk[k], tws[k], dih[k], ptws[k], rfy[k], tkf[k])
            crv_list.append(wp.nurbs_rep())

        ret = HorizontalStablizer()
        ret.surf = RuledSurf(crv_list[0], crv_list[1]) if len(crv_list) == 2 else Skinned(crv_list, 5, 3)
        ret.surf.pan(pan_dir)
        ret.tail = Coons(ret.surf.extract('U', 0), ret.surf.extract('U', 1), Line(ret.surf(0, 0), ret.surf(1, 0)), Line(ret.surf(0, 1), ret.surf(1, 1)))

        return ret


class VerticalStablizer(object):
    def __init__(self):
        self.surf = None
        self.tail = None

    @classmethod
    def from_section_param(cls, foil, zoff, cl, swp_bk, tws, dih, ptws, rfy, tkf, pan_dir, rot):
        """
        Construct the vertical stablizer from profile geometrical parameters.
        :param foil: Section airfoil list.
        :param zoff: Section 'z'-dim offset list when constructing from wing.
        :param cl: Section chord length list.
        :param swp_bk: Section sweep back angle list.
        :param tws:  Section twist angle list.
        :param dih: Section dihedral angle list.
        :param ptws: Section twist relative position list.
        :param rfy: Section 'y'-dim ref list when constructing from wing.
        :param tkf: Section thickness factor list.
        :param pan_dir: Panning direction before rotation.
        :param rot: Angle of rotation.
        :type rot: float
        :return: The vertical stablizer.
        :rtype: VerticalStablizer
        """

        '''Defensive check'''
        n = len(foil)
        if n < 2:
            raise ValueError("Insufficient input.")

        crv_list = []
        for k in range(n):
            wp = WingProfile.from_geom_param(foil[k], zoff[k], cl[k], swp_bk[k], tws[k], dih[k], ptws[k], rfy[k], tkf[k])
            crv_list.append(wp.nurbs_rep())

        ret = VerticalStablizer()
        ret.surf = RuledSurf(crv_list[0], crv_list[1]) if len(crv_list) == 2 else Skinned(crv_list, 5, 3)
        ret.surf.rotate((0, 0, 0), (-1, 0, 0), rot)
        ret.surf.pan(pan_dir)
        ret.tail = Coons(ret.surf.extract('U', 0), ret.surf.extract('U', 1), Line(ret.surf(0, 0), ret.surf(1, 0)), Line(ret.surf(0, 1), ret.surf(1, 1)))

        return ret


if __name__ == '__main__':
    unittest.main()
