import numpy as np
import math
from copy import deepcopy
from src.geom.utility import to_homogeneous, to_cartesian
from src.geom.curve import Arc, Line, ClampedNURBSCrv, ConicArc
from src.geom.surface import Coons

sqrt2 = math.sqrt(2)


class Nose(object):
    def __init__(self, nl, nfr, dh, h2, w2):
        """
        Simplified nose of a 'Wing-Tube' configuration aircraft.
        :param nl: Length of nose.
        :type nl: float
        :param nfr: Radius of the front hole.
        :type nfr: float
        :param dh:  Delta height between the center of the front section and the back section.
        :type dh: float
        :param h2: Half of the height of the back section.
        :type h2: float
        :param w2: Half of the width of the back section.
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

        arc3 = Arc.from_2pnt(p0, p3, 180, (1, 0, 0))
        arc5 = ConicArc(p6, t6, p8, t8, p7)
        for i in range(3):
            arc5.Pw[i][1] = i / 2 * p9[1] * arc5.Pw[i][-1]
        arc5.reset(arc5.U, arc5.Pw)

        '''Construct the ellipse by stretching the circle'''
        arc4 = Arc.from_2pnt(p2, p5, 180, (1, 0, 0))
        arc4.reset(arc4.U, np.copy(list(map(lambda u: to_homogeneous(to_cartesian(u) * (1, 1, b / a), u[-1]), arc4.Pw))))

        arc3_1, arc3_2 = ClampedNURBSCrv.split(arc3, [0.5])
        arc4_1, arc4_2 = ClampedNURBSCrv.split(arc4, [0.5])

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
        :param l1: Upper part of the front section curve.
        :type l1: Crv
        :param l2: Lower part of the front section curve.
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
        :param fu: Upper part of the front section.
        :type fu: Crv
        :param fd: Lower part of the front section.
        :type fd: Crv
        :param tl: Length of the tail.
        :type tl: float
        :param dh: Delta of the upper frame on back section.
        :type dh: float
        :param ttr: Radius of the hole on back section.
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

        arc8 = Arc.from_2pnt(p12, p15, 180, (1, 0, 0))
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
