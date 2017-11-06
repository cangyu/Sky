import unittest
import math
import os
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from abc import ABCMeta, abstractmethod
from math import sin, cos, tan, radians, fabs, atan2
from scipy.integrate import romberg
from scipy.interpolate import make_interp_spline
from scipy.optimize import root
from smooth import Laplace2D
from fluent import XF_MSH, BCType
from tfi import LinearTFI2D, LinearTFI3D
from misc import pnt_dist
from spacing import hyperbolic_tangent, uniform, single_exponential, double_exponential
from settings import AIRFOIL_DIR
from nurbs import Crv, ConicArc, to_homogeneous, to_cartesian
from nurbs import LocalCubicInterpolatedCrv, Line, point_inverse, GlobalInterpolatedCrv, Arc, Skinned, RuledSurf, Surf, Coons

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


class WingFrame(metaclass=ABCMeta):
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
    def y_front(self, u):
        """
        Calculate the Y-coordinate on leading edge.
        :param u: Relative position parameter.
        :type u: float
        :return: Y-coordinate on leading edge.
        :rtype: float
        """

        pass

    @abstractmethod
    def y_tail(self, u):
        """
        Calculate the Y-coordinate on leading edge.
        :param u: Relative position parameter.
        :type u: float
        :return: Y-coordinate on leading edge.
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
    def __str__(self):
        pass


class BWBFrame(WingFrame):
    def __init__(self, c_root, c_mid, c_tip, b_mid, b_tip, alpha_mid, alpha_tip):
        """
        BWB(Blended Wing Body)构型飞行器参数化描述
        :param c_root: 根部弦长
        :type c_root: float
        :param c_mid: 中间弦长
        :type c_mid: float
        :param c_tip: 翼尖弦长
        :type c_tip: float
        :param b_mid: 内翼宽度
        :type b_mid: float
        :param b_tip: 机翼半展长
        :type b_tip: float
        :param alpha_mid: 中段平均后掠角
        :type alpha_mid: float
        :param alpha_tip: 翼尖平均后掠角
        :type alpha_tip: float
        """

        '''Basic geometric parameters'''
        self.Cr = c_root
        self.Cm = c_mid
        self.Ct = c_tip
        self.Bm = b_mid
        self.Bt = b_tip
        self.Am = alpha_mid
        self.At = alpha_tip

        '''Calculate pivots on each curve'''
        front_pnt = np.zeros((3, 3))
        tail_pnt = np.zeros((3, 3))
        tail_pnt[0][0] = self.Cr
        front_pnt[1][0] = self.Bm * tan(radians(self.Am))
        front_pnt[1][2] = self.Bm
        tail_pnt[1][0] = front_pnt[1][0] + self.Cm
        tail_pnt[1][2] = self.Bm
        front_pnt[2][0] = front_pnt[1][0] + (self.Bt - self.Bm) * tan(radians(self.At))
        front_pnt[2][2] = self.Bt
        tail_pnt[2][0] = front_pnt[2][0] + self.Ct
        tail_pnt[2][2] = self.Bt

        '''Build interpolated functions'''
        self.u = np.array([0, self.Bm / self.Bt, 1.0])
        self.zw = np.array([0, self.Bm, self.Bt])
        self.xf = make_interp_spline(self.u, front_pnt[:, 0], 3, bc_type=([(1, 0)], [(2, 0)]))
        self.xt = make_interp_spline(self.u, tail_pnt[:, 0], 3, bc_type=([(1, 0)], [(2, 0)]))

    def x_front(self, u):
        return self.x_front(u)

    def x_tail(self, u):
        return self.xt(u)

    def y_front(self, u):
        return 0

    def y_tail(self, u):
        return 0

    def z(self, u):
        if u <= self.u[1]:
            return self.zw[1] * u / self.u[1]
        else:
            return self.zw[1] + (self.zw[2] - self.zw[1]) * (u - self.u[1]) / (self.u[2] - self.u[1])

    @property
    def span(self):
        return 2 * self.Bt

    @property
    def root_chord_len(self):
        return self.Cr

    @property
    def tip_chord_len(self):
        return self.Ct

    def __str__(self):
        a0 = self.area
        a1 = self.mean_aerodynamic_chord
        a2 = self.span
        a3 = self.root_chord_len
        a4 = self.tip_chord_len
        a5 = a3 / a4
        a6 = a2 ** 2 / a0

        ret = "Blended-Wing-Body Parametric Info:\n"
        ret += "General Info:\n"
        ret += "Area: {:.4f} m^2\n".format(a0)
        ret += "MAC: {:.3f} m\n".format(a1)
        ret += "Span: {:.3f} m\n".format(a2)
        ret += "Root Chord: {:.3f} m\n".format(a3)
        ret += "Tip Chord: {:.3f} m\n".format(a4)
        ret += "Taper Ratio: {:.3f}\n".format(a5)
        ret += "Aspect Ratio: {:.3f}\n".format(a6)
        ret += "Unique Info:\n"
        ret += "Inner Span {:.3f} m\n".format(self.Bm)
        ret += "Mid Chord: {:.3f} m\n".format(self.Cm)
        ret += "Inner wing sweep-back: {:.2f} (deg)\n".format(self.Am)
        ret += "Outer wing sweep-back: {:.2f} (deg)".format(self.At)

        return ret

    def show(self, section=None, n=1000):
        """
        Plot the leading edge, trailing edge, and profile dashes.
        :param section: Section distribution.
        :param n: Number of sampling points.
        :param n: int
        :return: None.
        """

        '''Show leading and trailing edge'''
        u_dist = np.linspace(0, 1.0, n)
        z = np.empty(n, float)
        xf = np.empty(n, float)
        xt = np.empty(n, float)
        for k in range(n):
            z[k] = self.z(u_dist[k])
            xf[k] = self.x_front(u_dist[k])
            xt[k] = self.x_tail(u_dist[k])

        plt.figure()
        plt.plot(z, xf, label='Leading Edge')
        plt.plot(z, xt, label='Trailing Edge')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.gca().set_aspect('equal')

        '''Show profiles'''
        if section is not None:
            for u in section:
                zp = self.z(u)
                plt.plot([zp, zp], [self.x_front(u), self.x_tail(u)], '--')

        plt.show()

    @classmethod
    def from_area_mac_span(cls, area2, mac, span2, c_root, c_tip, alpha_mid, alpha_tip):
        """
        Construct the planar frame with constant Area and MAC.
        :param area2: Half of the area of the wing.
        :type area2: float
        :param mac: Mean aerodynamic chord length of the wing.
        :type mac: float
        :param span2: Half of the width of the span.
        :type span2: float
        :param c_root: Length of the root chord.
        :type c_root: float
        :param c_tip: Length of the tip chord.
        :type c_tip: float
        :param alpha_mid: Averaged sweep back angle of the inner wing.
        :type alpha_mid: float
        :param alpha_tip: Averaged sweep back angle of the outer wing.
        :type alpha_tip: float
        :return: Constrained frame.
        :rtype: BWBFrame
        """

        '''Local constants'''
        tg1 = tan(radians(alpha_mid))
        tg2 = tan(radians(alpha_tip))

        def f(_x):
            """
            Constrained function.
            f0(_x) - area2 = 0
            f1(_x) - mac = 0
            :param _x: Variable vector.
            :return: LHS of the f(_x).
            """

            _bm = _x[0]
            _cm = _x[1]

            u = np.array([0, _bm / span2, 1.])
            fcx = np.array([0, _bm * tg1, _bm * tg1 + (span2 - _bm) * tg2])
            tcx = np.array([fcx[0] + c_root, fcx[1] + _cm, fcx[2] + c_tip])

            xf = make_interp_spline(u, fcx, 3, bc_type=([(1, 0)], [(2, 0)]))
            xt = make_interp_spline(u, tcx, 3, bc_type=([(1, 0)], [(2, 0)]))

            cur_s2 = span2 * romberg(lambda _u: xt(_u) - xf(_u), 0, 1)
            cur_mac = span2 * romberg(lambda _u: (xt(_u) - xf(_u)) ** 2, 0, 1) / cur_s2

            return cur_s2 - area2, cur_mac - mac

        '''Empirical guess'''
        pit = np.array([0.3 * span2, 0.5 * c_root])
        try:
            ans = root(f, pit)
            b_mid = ans.x[0]
            c_mid = ans.x[1]
        except:
            return None

        '''Build'''
        return BWBFrame(c_root, c_mid, c_tip, b_mid, span2, alpha_mid, alpha_tip)


class HWBFrame(WingFrame):
    def __init__(self, spn, cr, ct, fl, alpha, tl, beta, outer_taper=2.5):
        """
        Parametric wing Planform for HWB configuration.
        :param spn: Half of the span of the wing.
        :type spn: float
        :param cr: Length of the root-chord.
        :type cr: float
        :param ct: Length of the tip-chord.
        :type ct: float
        :param fl: Front segment length.
        :param alpha: Front segment sweep-back.
        :param tl: Tail segment length.
        :param beta: Tail segment sweep-back.
        :param outer_taper: Taper ratio of the outer wing.
        :type outer_taper: float
        """

        if len(fl) != 4:
            raise AssertionError("Invalid input of front segment length.")
        if len(alpha) != 4:
            raise AssertionError("Invalid input of front segment sweep-back.")
        if len(tl) != 4:
            raise AssertionError("Invalid input of tail segment length.")
        if len(beta) != 4:
            raise AssertionError("Invalid input of tail segment sweep-back.")

        '''Basic Wing Planform parameters'''
        self.Spn2 = spn
        self.Cr = cr
        self.Ct = ct
        self.OuterTaperRatio = outer_taper

        '''Front'''
        fl[-1] = spn - sum(fl[:-1])
        for l in fl:
            if l < 0:
                raise ValueError("Invalid parameter.")

        fp = np.zeros((5, 3))
        for i in range(4):
            fp[i + 1] = fp[i]
            fp[i + 1][0] += fl[i] * tan(alpha[i])
            fp[i + 1][2] += fl[i]
        ftv = np.array([[0, 0, 1],
                        [sin(alpha[1]), 0, cos(alpha[1])],
                        [sin(alpha[1]), 0, cos(alpha[1])],
                        [sin(alpha[3]), 0, cos(alpha[3])],
                        [sin(alpha[3]), 0, cos(alpha[3])]])

        '''Tail'''
        tl[-1] = fl[-1]
        tl[-2] = spn - (tl[0] + tl[1] + tl[-1])
        for l in tl:
            if l < 0:
                raise ValueError("Invalid parameter.")

        tp = np.zeros((5, 3))
        tp[0][0] = cr
        for i in range(2):
            tp[i + 1] = tp[i]
            tp[i + 1][0] += tl[i] * tan(beta[i])
            tp[i + 1][2] += tl[i]

        tp[-1] = fp[-1]
        tp[-1][0] += ct

        tp[-2] = fp[-2]
        tp[-2][0] += outer_taper * ct

        beta[-1] = atan2(tp[-1][0] - tp[-2][0], tl[-1])
        beta[-2] = atan2(tp[-2][0] - tp[-3][0], tl[-2])

        ttv = np.array([[0, 0, 1],
                        [sin(beta[1]), 0, cos(beta[1])],
                        [sin(beta[1]), 0, cos(beta[1])],
                        [sin(beta[3]), 0, cos(beta[3])],
                        [sin(beta[3]), 0, cos(beta[3])]])

        self.fl = np.copy(fl)
        self.alpha = np.copy(alpha)
        self.tl = np.copy(tl)
        self.beta = np.copy(beta)

        '''Interpolation'''
        self.front_crv = LocalCubicInterpolatedCrv(fp, ftv)
        self.tail_crv = LocalCubicInterpolatedCrv(tp, ttv)
        self.root_line = Line(fp[0], tp[0])
        self.tip_line = Line(fp[-1], tp[-1])
        self.planform_surf = Coons(self.root_line, self.tip_line, self.front_crv, self.tail_crv)

    def x_front(self, u):
        zp = self.z(u)
        lu = point_inverse(self.front_crv, zp, 2)
        return self.front_crv(lu)[0]

    def x_tail(self, u):
        zp = self.z(u)
        lu = point_inverse(self.tail_crv, zp, 2)
        return self.tail_crv(lu)[0]

    def y_front(self, u):
        return 0

    def y_tail(self, u):
        return 0

    def z(self, u):
        return u * self.Spn2

    @property
    def span(self):
        return 2 * self.Spn2

    @property
    def root_chord_len(self):
        return self.Cr

    @property
    def tip_chord_len(self):
        return self.Ct

    def __str__(self):
        a0 = self.area
        a1 = self.mean_aerodynamic_chord
        a2 = self.span
        a3 = self.root_chord_len
        a4 = self.tip_chord_len
        a5 = a3 / a4
        a6 = a2 ** 2 / a0

        ret = "HWB Wing Planform Info:\n"
        ret += "General:\n"
        ret += "Area: {:.4f} m^2\n".format(a0)
        ret += "MAC: {:.3f} m\n".format(a1)
        ret += "Span: {:.3f} m\n".format(a2)
        ret += "Root: {:.3f} m\n".format(a3)
        ret += "Tip: {:.3f} m\n".format(a4)
        ret += "Taper Ratio: {:.3f}\n".format(a5)
        ret += "AR: {:.3f}\n".format(a6)
        ret += "Unique:\n"
        ret += "Outer Span: {} m\n".format(self.fl[-1])
        ret += "Outer AR: {}\n".format(2 * self.fl[-1] / (self.Ct * (1 + self.OuterTaperRatio) / 2))

        return ret

    def show(self, dist, n=1000):
        u_dist = np.linspace(0, 1.0, n)
        zf = np.empty(n, float)
        zt = np.empty(n, float)
        xf = np.empty(n, float)
        xt = np.empty(n, float)
        for k in range(n):
            fp = self.front_crv(u_dist[k])
            tp = self.tail_crv(u_dist[k])
            xf[k] = fp[0]
            zf[k] = fp[2]
            xt[k] = tp[0]
            zt[k] = tp[2]

        spn = self.front_crv.end[2]
        cdst = np.copy(dist) * spn
        fu = np.copy(cdst)
        tu = np.copy(cdst)
        for k, u in enumerate(fu):
            fu[k] = point_inverse(self.front_crv, u, 2)
        for k, u in enumerate(tu):
            tu[k] = point_inverse(self.tail_crv, u, 2)

        plt.figure()
        plt.plot(zf, xf, label='Leading Edge')
        plt.plot(zt, xt, label='Trailing Edge')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.gca().set_aspect('equal')

        for k, u in enumerate(cdst):
            tfx = self.front_crv(fu[k])[0]
            ttx = self.tail_crv(tu[k])[0]
            plt.plot([u, u], [tfx, ttx], '--')
            plt.text(u, (tfx + ttx) / 2, str(k))

        plt.show()


class Airfoil(object):
    AIRFOIL_LIST = []

    def __init__(self):
        """
        2D Airfoil, with chord length equals to 1.
        """

        self.name = None
        self.pts = None

    def __repr__(self):
        ret = "Airfoil: {}\nNbrOfPnt: {}\nCoordinates(from tail-up to tail-down):\n".format(self.name, len(self.pts))
        crd = []
        for k, pnt in enumerate(self.pts):
            crd.append("{}: ({}, {}, {})".format(k + 1, pnt[0], pnt[1], pnt[2]))

        ret += "\n".join(crd)
        return ret

    def __str__(self):
        return self.name

    def to_blunt(self):
        r = self.pts[-2][0]
        for k, p in enumerate(self.pts):
            self.pts[k][0] = p[0] / r
        self.pts = self.pts[1:-1]

    @classmethod
    def update_airfoil_list(cls):
        for f in os.listdir(AIRFOIL_DIR):
            base, ext = os.path.splitext(f)
            if ext == '.dat':
                cls.AIRFOIL_LIST.append(base)

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
    def tail(self):
        return (self.pts[0] + self.pts[-1]) / 2

    @property
    def tail_up(self):
        return self.pts[0]

    @property
    def tail_down(self):
        return self.pts[-1]

    @property
    def chord_len(self):
        return pnt_dist(self.front, self.tail)

    @property
    def is_blunt(self):
        return not math.isclose(norm(self.pts[0] - self.pts[-1], np.inf), 0.0)

    @property
    def pnt_num(self):
        return len(self.pts)

    def show(self):
        (px, py, pz) = zip(*self.pts)
        plt.plot(px, py)
        plt.gca().set_aspect('equal')
        plt.show()

    def plot_curvature(self, n=2000):
        """
        Plot the curvature along the airfoil from tail-up to tail_down, anti-clockwise.
        By default, the interpolated curve is cubic
        :param n: Number of sampling points.
        :type n: int
        :return: None
        """

        '''Build curve'''
        crv = self.nurbs_rep()

        '''Calculate curvature'''
        ul = np.linspace(crv.U[0], crv.U[-1], n)
        kappa = list(map(lambda u: crv.curvature(u), ul))

        '''Plot'''
        plt.figure()
        plt.plot(ul, kappa)
        plt.xlabel('Parameter along curve')
        plt.ylabel('Curvature')
        plt.title(self.name)
        plt.show()

    def nurbs_rep(self, p=3, method='centripetal'):
        """
        NURBS Representation
        :param p: 插值次数
        :type p: int
        :param method: 插值参数化方法
        :type method: str
        :return: 翼型的NURBS全局插值曲线
        :rtype: GlobalInterpolatedCrv
        """

        return GlobalInterpolatedCrv(self.pts, p, method)

    @classmethod
    def from_database(cls, fn):
        """
        Construct airfoil from local database.
        :param fn: Name of the airfoil. (Capital case).
        :type fn: str
        :return: Airfoil Object.
        :rtype: Airfoil
        """

        af = cls()
        af.name = fn
        af.pts = cls.read_pts(fn)
        return af

    def gen_grid(self):
        crv = self.nurbs_rep()

        '''Vertical boundary'''
        calc_dir = lambda start, end: (end[0] - start[0], end[1] - start[1])
        calc_len = lambda dir: math.pow(math.pow(dir[0], 2) + math.pow(dir[1], 2), 0.5)
        inner_product = lambda d1, d2: d1[0] * d2[0] + d1[1] * d2[1]

        dir1 = calc_dir(af.pts[1], af.pts[0])
        len1 = calc_len(dir1)
        dir2 = calc_dir(af.pts[-2], af.pts[-1])
        len2 = calc_len(dir2)

        rotate = lambda dir, theta: (dir[0] * math.cos(math.radians(theta)) - dir[1] * math.sin(math.radians(theta)),
                                     dir[0] * math.sin(math.radians(theta)) + dir[1] * math.cos(math.radians(theta)))
        dir1 = rotate(dir1, 90)
        dir2 = rotate(dir2, 270)
        dir3 = (dir1[0] / len1, dir1[1] / len1)
        dir4 = (dir2[0] / len2, dir2[1] / len2)

        c2 = lambda v: np.array([af.pts[0][0] + v * R * dir3[0], af.pts[0][1] + v * R * dir3[1], 0])
        r = calc_len(c2(1.0))

        dir5 = (af.pts[-1][0], af.pts[-1][1])
        l4 = calc_len(dir5)
        theta = math.pi - math.acos(inner_product(dir4, dir5) / l4)
        alpha = math.asin(l4 / r * math.sin(theta))
        beta = math.pi - theta - alpha
        b = r / math.sin(theta) * math.sin(beta)

        c4 = lambda v: np.array([af.pts[-1][0] + v * b * dir4[0], af.pts[-1][1] + v * b * dir4[1], 0])

        '''Farfield boundary'''
        sp = c2(1.0)
        ep = c4(1.0)
        sa = math.atan2(sp[1], sp[0])
        ea = math.atan2(ep[1], ep[0])
        if ea < 0:
            ea += math.pi * 2

        c3 = lambda u: np.array([r * math.cos((1 - u) * sa + u * ea), r * math.sin((1 - u) * sa + u * ea), 0])

        '''Tail'''
        c5 = lambda u: np.array([(1 - u) * af.pts[-1][0] + u * af.pts[0][0], (1 - u) * af.pts[-1][1] + u * af.pts[0][1], 0])
        ea2 = ea - math.pi * 2 if ea > 0 else ea
        c6 = lambda u: np.array([r * math.cos((1 - u) * ea2 + u * sa), r * math.sin((1 - u) * ea2 + u * sa), 0])

        # return c1, c2, c3, c4, c5, c6

        def write_uniform_airfoil(foil, L, R, U1, U2, V, fn="", delta_zeta=1.0, delta_eta=1.0):
            u1_list = np.linspace(0, 1.0, U1 + 1)
            u2_list = np.linspace(0, 1.0, U2 + 1)
            v_list = np.linspace(0, 1.0, V + 1)
            c1, c2, c3, c4, c5, c6 = airfoil(foil, L, R)

            grid1 = LinearTFI2D(c1, c2, c3, c4)
            ppu, ppv = np.meshgrid(u1_list, v_list, indexing='ij')
            grid1.calc_grid(ppu, ppv)
            grid1 = Laplace2D(grid1.grid)
            grid1.smooth()

            grid2 = LinearTFI2D(c2, c5, c4, c6)
            ppu, ppv = np.meshgrid(u2_list, v_list, indexing='ij')
            grid2.calc_grid(ppu, ppv)

            if fn == "":
                fn += foil
                fn += "_{}_{}_{}_{}_{}_{}_{}_Multi.xyz".format(L, R, U1, U2, V, delta_zeta, delta_eta)

            msh = Plot3D()
            msh.add_block(PLOT3D_Block.build_from_2d(grid1.grid))
            msh.add_block(PLOT3D_Block.build_from_2d(grid2.grid))
            msh.write(fn)


class WingProfile(Airfoil):
    SEC_INTRINSIC_PARAM = ['Airfoil', 'Z(m)', 'X_front(m)', 'Y_front(m)', 'X_tail(m)', 'Y_tail(m)', 'Thickness Ratio']
    SEC_GEOM_PARAM = ['Airfoil', 'Z(m)', 'Length(m)', 'SweepBack(deg)', 'Twist(deg)', 'Dihedral(deg)', 'TwistPos', 'Thickness Ratio']

    def __init__(self, foil, ends, thickness_factor=1.0):
        """
        3D profile at certain position.
        :param foil: Airfoil name.
        :type foil: str
        :param ends: Starting and ending points of the section.
        :param thickness_factor: Vertical stretching factor.
        :type thickness_factor: float
        """

        super(WingProfile, self).__init__()

        '''Inspect endings'''
        self.ending = np.copy(ends)
        if not equal(ends[0][2], ends[1][2]):
            raise AssertionError("Invalid ending coordinates in Z direction!")

        cl = self.chord_len
        if equal(cl, 0.0):
            raise ZeroDivisionError("Invalid ending coordinates in XY direction!")

        rotation = complex((ends[1][0] - ends[0][0]) / cl, (ends[1][1] - ends[0][1]) / cl)

        '''Build section'''
        self.name = foil
        tfoil = Airfoil().from_database(foil)
        if not tfoil.is_blunt:
            tfoil.to_blunt()
        self.pts = np.copy(tfoil.pts)
        n = self.pnt_num
        for i in range(n):
            '''Stretch, Z offset and Thickness'''
            self.pts[i][0] *= cl
            self.pts[i][1] *= (cl * thickness_factor)
            self.pts[i][2] = ends[0][2]

            '''Rotate around ends[0]'''
            origin_vector = complex(self.pts[i][0], self.pts[i][1])
            origin_vector *= rotation
            self.pts[i][0] = origin_vector.real
            self.pts[i][1] = origin_vector.imag

            '''Move to ends[0]'''
            self.pts[i][0] += ends[0][0]
            self.pts[i][1] += ends[0][1]

    @property
    def front(self):
        return self.ending[0]

    @property
    def tail(self):
        return self.ending[-1]

    @property
    def z_offset(self):
        return self.ending[0][-1]

    @classmethod
    def from_geom_param(cls, foil, z_offset, length, sweep_back, twist, dihedral, twist_pos=0.25, y_ref=0, thickness_factor=1.0):
        """
        从几何描述参数构建机翼剖面
        :param foil: 翼型名称
        :type foil: str
        :param z_offset: Z方向偏移量
        :type z_offset: float
        :param length: 剖面长度
        :type length: float
        :param sweep_back: 后掠角
        :type sweep_back: float
        :param twist: 相对翼根弦线的扭转角(默认在1/4弦长处扭转)
        :type twist: float
        :param dihedral: 相对翼根的上反角
        :type dihedral: float
        :param twist_pos: 扭转中心
        :type twist_pos: float
        :param y_ref: 翼根处Y方向基准坐标
        :type y_ref: float
        :param thickness_factor: 纵向厚度拉伸系数
        :type thickness_factor: float
        :return: 机翼剖面
        :rtype: WingProfile
        """

        x_offset = z_offset * math.tan(math.radians(sweep_back))
        y_offset = y_ref + z_offset * math.tan(math.radians(dihedral))
        front = np.array([x_offset, y_offset, z_offset], float)
        tail = np.array([x_offset + length, y_offset, z_offset], float)

        center = (1 - twist_pos) * front + twist_pos * tail
        theta = math.radians(-twist)
        rot = complex(math.cos(theta), math.sin(theta))
        d1 = front - center
        d2 = tail - center
        c1 = complex(d1[0], d1[1]) * rot
        c2 = complex(d2[0], d2[1]) * rot
        front[0] = center[0] + c1.real
        front[1] = center[1] + c1.imag
        tail[0] = center[0] + c2.real
        tail[1] = center[1] + c2.imag
        ending = np.empty((2, 3), float)
        ending[0] = front
        ending[1] = tail

        return cls(foil, ending, thickness_factor)


class Wing(object):
    def __init__(self, section_list):
        """
        从剖面序列构造机翼
        :param section_list: 机翼剖面序列
        """

        self.section = deepcopy(section_list)

    def __repr__(self):
        return "Wing with {} sections".format(self.section_num)

    @property
    def section_num(self):
        return len(self.section)

    @property
    def root(self):
        return self.section[0].nurbs_rep()

    @property
    def tip(self):
        return self.section[-1].nurbs_rep()

    def front(self, q=3, method='chord'):
        n = self.section_num
        front_pts = np.zeros((n, 3))
        for i in range(n):
            front_pts[i] = self.section[i].front

        return GlobalInterpolatedCrv(front_pts, q, method)

    def surf(self, p=5, q=3):
        """
        构建机翼轮廓曲线、曲面
        :param p: U方向次数
        :type p: int
        :param q: V方向次数
        :type q: int
        :return: 机翼蒙皮曲面
        :rtype: Skinned
        """

        profile_list = []
        for elem in self.section:
            profile_list.append(elem.nurbs_rep(p))
        return Skinned(profile_list, p, q)

    @property
    def tail_up(self):
        sk = self.surf()
        return sk.extract('U', 0)

    @property
    def tail_down(self):
        sk = self.surf()
        return sk.extract('U', 1)

    def iges_model(self, fn, p=3, q=3, mirror=True):
        """
        生成机翼相应的IGES模型
        :param fn: 文件名
        :type fn: str
        :param p: U方向次数
        :type p: int
        :param q: V方向次数
        :type q: int
        :param mirror: 是否生成对称部分
        :type mirror: bool
        :return: 可用于后续生成IGS文件的IGES_Model对象
        :rtype: IGES_Model
        """

        wing_model = IGES_Model(fn)

        '''前缘曲线'''
        wing_model.add_entity(self.front(q).to_iges())

        '''剖面'''
        for elem in self.section:
            wing_model.add_entity(elem.nurbs_rep(p).to_iges())

        '''蒙皮'''
        sk = self.surf(p, q)
        wing_model.add_entity(sk.to_iges())

        '''镜像'''
        if mirror:
            msk = deepcopy(sk)
            for i in range(msk.n + 1):
                for j in range(msk.m + 1):
                    msk.Pw[i][j][2] *= -1

            wing_model.add_entity(msk.to_iges())

        return wing_model

    @classmethod
    def from_intrinsic_desc(cls, airfoil, thickness, z, xf, yf, xt, yt):
        n = len(airfoil)
        section_list = []
        for k in range(n):
            ends = np.empty((2, 3), float)
            ends[0][2] = ends[1][2] = z[k]
            ends[0][0] = xf[k]
            ends[0][1] = yf[k]
            ends[1][0] = xt[k]
            ends[1][1] = yt[k]
            section_list.append(WingProfile(airfoil[k], ends, thickness[k]))

        return cls(section_list)

    @classmethod
    def from_geom_desc(cls, airfoil, length, thickness, z, sweep, twist, twist_pos, dihedral, y_ref):
        n = len(airfoil)
        section_list = []
        for k in range(n):
            section_list.append(WingProfile.from_geom_param(airfoil[k], z[k], length[k], sweep[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness[k]))

        return cls(section_list)

    @classmethod
    def from_frame(cls, airfoil, thickness, u, frm):
        """
        根据给定的参数化模型生成机翼
        :param airfoil: 剖面翼型序列
        :param thickness: 剖面厚度拉伸系数
        :param u: 剖面位置分布参数
        :param frm: 参数化模型
        :type frm: WingFrame
        :return:
        """

        z = list(map(frm.z, u))
        xf = list(map(frm.x_front, u))
        yf = list(map(frm.y_front, u))
        xt = list(map(frm.x_tail, u))
        yt = list(map(frm.y_tail, u))
        return cls.from_intrinsic_desc(airfoil, thickness, z, xf, yf, xt, yt)

    def gen_grid(self, fn, enn=50, inner_brk=(0.44, 0.56), outer_brk=(0.3, 0.72)):
        """
        Generate the multi-block grid for a simplified wing.
        :return: None
        """

        '''Grid parameters'''
        la = self.section[0].chord_len
        lt = 30 * la
        r = 10 * la
        inner_spn = self.section[-1].z_offset
        outer_spn = 20 * inner_spn

        wsf = self.surf()
        crv_root = wsf.extract('V', 0)
        crv_tip = wsf.extract('V', 1)
        far = WingProfile.from_geom_param(self.section[-1].name, inner_spn + outer_spn, la, 0, 0, 0, thickness_factor=3)
        crv_far = far.nurbs_rep()
        fsf = RuledSurf(crv_tip, crv_far)

        brk_root = brk_tip = brk_far = inner_brk
        obrk_root = obrk_tip = obrk_far = outer_brk

        '''Points, lines, curves, surfs'''
        p = np.zeros((36, 3))
        p[2] = wsf(0, 0)
        p[4] = wsf(1, 0)
        p[0] = p[2]
        p[0][1] += r
        p[6] = p[4]
        p[6][1] -= r
        p[1] = p[0]
        p[1][0] += lt
        p[3] = p[2]
        p[3][0] += lt
        p[5] = p[4]
        p[5][0] += lt
        p[7] = p[6]
        p[7][0] += lt
        p[9] = p[1]
        p[9][2] += inner_spn
        p[11] = p[3]
        p[11][2] += inner_spn
        p[13] = p[5]
        p[13][2] += inner_spn
        p[15] = p[7]
        p[15][2] += inner_spn
        p[8] = p[0]
        p[8][2] += inner_spn
        p[10] = wsf(0, 1)
        p[12] = wsf(1, 1)
        p[14] = p[6]
        p[14][2] += inner_spn
        p[16] = p[8]
        p[16][2] += outer_spn
        p[17] = p[9]
        p[17][2] += outer_spn
        p[22] = p[14]
        p[22][2] += outer_spn
        p[23] = p[15]
        p[23][2] += outer_spn
        p[18] = crv_far.start
        p[20] = crv_far.end
        p[19] = p[18]
        p[19][0] += lt
        p[21] = p[20]
        p[21][0] += lt
        p[24] = crv_root(brk_root[0])
        p[25] = crv_root(brk_root[1])
        p[26] = crv_tip(brk_tip[0])
        p[27] = crv_tip(brk_tip[1])
        p[28] = crv_far(brk_far[0])
        p[29] = crv_far(brk_far[1])
        outer_root = Arc.from_2pnt(p[0], p[6], 180, (0, 0, 1))
        outer_tip = Arc.from_2pnt(p[8], p[14], 180, (0, 0, 1))
        outer_far = Arc.from_2pnt(p[16], p[22], 180, (0, 0, 1))
        p[30] = outer_root(obrk_root[0])
        p[31] = outer_root(obrk_root[1])
        p[32] = outer_tip(obrk_tip[0])
        p[33] = outer_tip(obrk_tip[1])
        p[34] = outer_far(obrk_far[0])
        p[35] = outer_far(obrk_far[1])

        l = [Line(p[0], p[1]),  # 0
             Line(p[2], p[3]),  # 1
             Line(p[4], p[5]),  # 2
             Line(p[6], p[7]),  # 3
             Line(p[8], p[9]),  # 4
             Line(p[10], p[11]),  # 5
             Line(p[12], p[13]),  # 6
             Line(p[14], p[15]),  # 7
             Line(p[16], p[17]),  # 8
             Line(p[18], p[19]),  # 9
             Line(p[20], p[21]),  # 10
             Line(p[22], p[23]),  # 11
             Line(p[1], p[9]),  # 12
             Line(p[3], p[11]),  # 13
             Line(p[5], p[13]),  # 14
             Line(p[7], p[15]),  # 15
             Line(p[9], p[17]),  # 16
             Line(p[11], p[19]),  # 17
             Line(p[13], p[21]),  # 18
             Line(p[15], p[23]),  # 19
             Line(p[0], p[8]),  # 20
             Line(p[6], p[14]),  # 21
             Line(p[8], p[16]),  # 22
             Line(p[10], p[18]),  # 23
             Line(p[12], p[20]),  # 24
             Line(p[14], p[22]),  # 25
             Line(p[3], p[1]),  # 26
             Line(p[5], p[3]),  # 27
             Line(p[5], p[7]),  # 28
             Line(p[11], p[9]),  # 29
             Line(p[13], p[11]),  # 30
             Line(p[13], p[15]),  # 31
             Line(p[19], p[17]),  # 32
             Line(p[21], p[19]),  # 33
             Line(p[21], p[23]),  # 34
             Line(p[2], p[0]),  # 35
             Line(p[4], p[6]),  # 36
             Line(p[10], p[8]),  # 37
             Line(p[12], p[14]),  # 38
             Line(p[4], p[2]),  # 39
             Line(p[12], p[10]),  # 40
             Line(p[18], p[16]),  # 41
             Line(p[20], p[18]),  # 42
             Line(p[20], p[22]),  # 43
             Line(p[24], p[30]),  # 44
             Line(p[25], p[31]),  # 45
             Line(p[26], p[32]),  # 46
             Line(p[27], p[33]),  # 47
             Line(p[28], p[34]),  # 48
             Line(p[29], p[35]),  # 49
             Line(p[26], p[28]),  # 50
             Line(p[27], p[29]),  # 51
             Line(p[30], p[32]),  # 52
             Line(p[31], p[33]),  # 53
             Line(p[32], p[34]),  # 54
             Line(p[33], p[35])  # 55
             ]

        c = [wsf.extract('U', 0), wsf.extract('U', 1)]
        c2, c3, c4 = ClampedNURBSCrv.split(outer_root, obrk_root)
        c5, c6, c7 = ClampedNURBSCrv.split(outer_tip, obrk_tip)
        c8, c9, c10 = ClampedNURBSCrv.split(outer_far, obrk_far)
        c11, c12, c13 = ClampedNURBSCrv.split(crv_root, brk_root)
        c14, c15, c16 = ClampedNURBSCrv.split(crv_tip, brk_tip)
        c17, c18, c19 = ClampedNURBSCrv.split(crv_far, brk_far)
        c20 = wsf.extract('U', brk_root[0])
        c21 = wsf.extract('U', brk_root[1])
        c22 = fsf.extract('U', brk_tip[0])
        c23 = fsf.extract('U', brk_tip[1])
        c3.reverse()
        c4.reverse()
        c6.reverse()
        c7.reverse()
        c9.reverse()
        c10.reverse()
        c12.reverse()
        c13.reverse()
        c15.reverse()
        c16.reverse()
        c18.reverse()
        c19.reverse()
        c.append(c2)
        c.append(c3)
        c.append(c4)
        c.append(c5)
        c.append(c6)
        c.append(c7)
        c.append(c8)
        c.append(c9)
        c.append(c10)
        c.append(c11)
        c.append(c12)
        c.append(c13)
        c.append(c14)
        c.append(c15)
        c.append(c16)
        c.append(c17)
        c.append(c18)
        c.append(c19)
        c.append(c20)
        c.append(c21)
        c.append(c22)
        c.append(c23)

        ts1 = Surf.split(wsf, brk_root, [])
        s0 = ts1[0][0]
        s1 = ts1[1][0]
        s2 = ts1[2][0]
        ts2 = Surf.split(fsf, brk_tip, [])
        s3 = ts2[0][0]
        s4 = ts2[1][0]
        s5 = ts2[2][0]
        s = [s0, s1, s2, s3, s4, s5]

        '''Node number and distribution'''
        n = np.full(8, enn, int)
        u0 = hyperbolic_tangent(n[0], 8)
        u1 = double_exponential(n[1], 0.5, 1.5, 0.5)
        u2 = uniform(n[2])
        u3 = single_exponential(n[3], 5)
        u4 = hyperbolic_tangent(n[4], 5)
        u5 = double_exponential(n[5], 0.5, 1.2, 0.5)
        u6 = double_exponential(n[6], 0.5, 1.5, 0.5)
        u7 = uniform(n[7])
        knot_dist = [u0, u1, u2, u3, u4, u5, u6, u7]

        '''Construct blocks'''
        blk_list = []
        blk_param_list = []

        b0_tfi_grid = LinearTFI3D.from_edges(l[1], l[26], l[0], l[35], l[5], l[29], l[4], l[37], c[0], l[13], l[12], l[20])
        blk_list.append(b0_tfi_grid)
        blk_param_list.append([knot_dist[3], knot_dist[0], knot_dist[2]])

        b1_tfi_grid = LinearTFI3D.from_edges(l[2], l[27], l[1], l[39], l[6], l[30], l[5], l[40], c[1], l[14], l[13], c[0])
        blk_list.append(b1_tfi_grid)
        blk_param_list.append([knot_dist[3], knot_dist[7], knot_dist[2]])

        b2_tfi_grid = LinearTFI3D.from_edges(l[36], l[3], l[28], l[2], l[38], l[7], l[31], l[6], c[1], l[21], l[15], l[14])
        blk_list.append(b2_tfi_grid)
        blk_param_list.append([knot_dist[0], knot_dist[3], knot_dist[2]])

        b3_tfi_grid = LinearTFI3D.from_edges(l[5], l[29], l[4], l[37], l[9], l[32], l[8], l[41], l[23], l[17], l[16], l[22])
        blk_list.append(b3_tfi_grid)
        blk_param_list.append([knot_dist[3], knot_dist[0], knot_dist[4]])

        b4_tfi_grid = LinearTFI3D.from_edges(l[6], l[30], l[5], l[40], l[10], l[33], l[9], l[42], l[24], l[18], l[17], l[23])
        blk_list.append(b4_tfi_grid)
        blk_param_list.append([knot_dist[3], knot_dist[7], knot_dist[4]])

        b5_tfi_grid = LinearTFI3D.from_edges(l[38], l[7], l[31], l[6], l[43], l[11], l[34], l[10], l[24], l[25], l[19], l[18])
        blk_list.append(b5_tfi_grid)
        blk_param_list.append([knot_dist[0], knot_dist[3], knot_dist[4]])

        b6_s1 = deepcopy(s[0])
        b6_s2 = LinearTFI2D(c[2], l[20], c[5], l[52])
        b6_s3 = LinearTFI2D(c[0], l[35], l[20], l[37])
        b6_s4 = LinearTFI2D(c[20], l[44], l[52], l[46])
        b6_s5 = LinearTFI2D(l[35], c[11], l[44], c[2])
        b6_s6 = LinearTFI2D(l[37], c[14], l[46], c[5])

        b6_tfi_grid = LinearTFI3D(lambda v, w: b6_s1(v, w),
                                  lambda v, w: b6_s2(v, w),
                                  lambda w, u: b6_s3(w, u),
                                  lambda w, u: b6_s4(w, u),
                                  lambda u, v: b6_s5(u, v),
                                  lambda u, v: b6_s6(u, v))

        blk_list.append(b6_tfi_grid)
        blk_param_list.append([knot_dist[0], knot_dist[1], knot_dist[2]])

        b7_s1 = LinearTFI2D(l[45], c[21], l[47], l[53])
        b7_s2 = LinearTFI2D(l[44], c[20], l[46], l[52])
        b7_s3 = deepcopy(s[1])
        b7_s3.reverse('U')
        b7_s3.swap()
        b7_s4 = LinearTFI2D(l[53], c[3], l[52], c[6])
        b7_s5 = LinearTFI2D(c[12], l[45], c[3], l[44])
        b7_s6 = LinearTFI2D(c[15], l[47], c[6], l[46])

        b7_tfi_grid = LinearTFI3D(lambda v, w: b7_s1(v, w),
                                  lambda v, w: b7_s2(v, w),
                                  lambda w, u: b7_s3(w, u),
                                  lambda w, u: b7_s4(w, u),
                                  lambda u, v: b7_s5(u, v),
                                  lambda u, v: b7_s6(u, v))

        blk_list.append(b7_tfi_grid)
        blk_param_list.append([knot_dist[5], knot_dist[0], knot_dist[2]])

        b8_s1 = LinearTFI2D(l[36], c[1], l[38], l[21])
        b8_s2 = LinearTFI2D(l[45], c[21], l[47], l[53])
        b8_s3 = deepcopy(s[2])
        b8_s3.reverse('U')
        b8_s3.swap()
        b8_s4 = LinearTFI2D(l[21], c[4], l[53], c[7])
        b8_s5 = LinearTFI2D(c[13], l[36], c[4], l[45])
        b8_s6 = LinearTFI2D(c[16], l[38], c[7], l[47])

        b8_tfi_grid = LinearTFI3D(lambda v, w: b8_s1(v, w),
                                  lambda v, w: b8_s2(v, w),
                                  lambda w, u: b8_s3(w, u),
                                  lambda w, u: b8_s4(w, u),
                                  lambda u, v: b8_s5(u, v),
                                  lambda u, v: b8_s6(u, v))

        blk_list.append(b8_tfi_grid)
        blk_param_list.append([knot_dist[6], knot_dist[0], knot_dist[2]])

        b9_s1 = deepcopy(s[3])
        b9_s2 = LinearTFI2D(c[5], l[22], c[8], l[54])
        b9_s3 = LinearTFI2D(l[23], l[37], l[22], l[41])
        b9_s4 = LinearTFI2D(c[22], l[46], l[54], l[48])
        b9_s5 = LinearTFI2D(l[37], c[14], l[46], c[5])
        b9_s6 = LinearTFI2D(l[41], c[17], l[48], c[8])

        b9_tfi_grid = LinearTFI3D(lambda v, w: b9_s1(v, w),
                                  lambda v, w: b9_s2(v, w),
                                  lambda w, u: b9_s3(w, u),
                                  lambda w, u: b9_s4(w, u),
                                  lambda u, v: b9_s5(u, v),
                                  lambda u, v: b9_s6(u, v))

        blk_list.append(b9_tfi_grid)
        blk_param_list.append([knot_dist[0], knot_dist[1], knot_dist[4]])

        b10_s1 = LinearTFI2D(l[47], c[23], l[49], l[55])
        b10_s2 = LinearTFI2D(l[46], c[22], l[48], l[54])
        b10_s3 = deepcopy(s[4])
        b10_s3.reverse('U')
        b10_s3.swap()
        b10_s4 = LinearTFI2D(l[55], c[6], l[54], c[9])
        b10_s5 = LinearTFI2D(c[15], l[47], c[6], l[46])
        b10_s6 = LinearTFI2D(c[18], l[49], c[9], l[48])

        b10_tfi_grid = LinearTFI3D(lambda v, w: b10_s1(v, w),
                                   lambda v, w: b10_s2(v, w),
                                   lambda w, u: b10_s3(w, u),
                                   lambda w, u: b10_s4(w, u),
                                   lambda u, v: b10_s5(u, v),
                                   lambda u, v: b10_s6(u, v))

        blk_list.append(b10_tfi_grid)
        blk_param_list.append([knot_dist[5], knot_dist[0], knot_dist[4]])

        b11_s1 = LinearTFI2D(l[38], l[24], l[43], l[25])
        b11_s2 = LinearTFI2D(l[47], c[23], l[49], l[55])
        b11_s3 = deepcopy(s[5])
        b11_s3.reverse('U')
        b11_s3.swap()
        b11_s4 = LinearTFI2D(l[25], c[7], l[55], c[10])
        b11_s5 = LinearTFI2D(c[16], l[38], c[7], l[47])
        b11_s6 = LinearTFI2D(c[19], l[43], c[10], l[49])

        b11_tfi_grid = LinearTFI3D(lambda v, w: b11_s1(v, w),
                                   lambda v, w: b11_s2(v, w),
                                   lambda w, u: b11_s3(w, u),
                                   lambda w, u: b11_s4(w, u),
                                   lambda u, v: b11_s5(u, v),
                                   lambda u, v: b11_s6(u, v))

        blk_list.append(b11_tfi_grid)
        blk_param_list.append([knot_dist[6], knot_dist[0], knot_dist[4]])

        b12_s1 = deepcopy(s[5])
        b12_s1.reverse('U')
        b12_s2 = deepcopy(s[3])
        b12_s3 = LinearTFI2D(l[24], l[40], l[23], l[42])
        b12_s4 = deepcopy(s[4])
        b12_s4.reverse('U')
        b12_s4.swap()
        b12_s5 = LinearTFI2D(l[40], c[16], c[15], c[14])
        b12_s6 = LinearTFI2D(l[42], c[19], c[18], c[17])

        b12_tfi_grid = LinearTFI3D(lambda v, w: b12_s1(v, w),
                                   lambda v, w: b12_s2(v, w),
                                   lambda w, u: b12_s3(w, u),
                                   lambda w, u: b12_s4(w, u),
                                   lambda u, v: b12_s5(u, v),
                                   lambda u, v: b12_s6(u, v))

        blk_list.append(b12_tfi_grid)
        blk_param_list.append([knot_dist[7], knot_dist[6], knot_dist[4]])

        def report(msg):
            print('Process {} : {}'.format(os.getpid(), msg))

        report('Calculating grid ...')
        for i in range(len(blk_list)):
            _u, _v, _w = blk_param_list[i]
            blk_list[i].calc_grid(_u, _v, _w)
            report('Block {} calculation done!'.format(i))

        '''网格, 边界条件, 邻接关系'''
        blk = [b0_tfi_grid.get_grid(),
               b1_tfi_grid.get_grid(),
               b2_tfi_grid.get_grid(),
               b3_tfi_grid.get_grid(),
               b4_tfi_grid.get_grid(),
               b5_tfi_grid.get_grid(),
               b6_tfi_grid.get_grid(),
               b7_tfi_grid.get_grid(),
               b8_tfi_grid.get_grid(),
               b9_tfi_grid.get_grid(),
               b10_tfi_grid.get_grid(),
               b11_tfi_grid.get_grid(),
               b12_tfi_grid.get_grid()]

        bc = [(BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b0
              (BCType.Wall, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Symmetry, BCType.Interior),  # b1
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b2
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b3
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField),  # b4
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b5
              (BCType.Wall, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Symmetry, BCType.Interior),  # b6
              (BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b7
              (BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField, BCType.Symmetry, BCType.Interior),  # b8
              (BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField),  # b9
              (BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b10
              (BCType.Interior, BCType.Interior, BCType.Interior, BCType.PressureFarField, BCType.Interior, BCType.PressureFarField),  # b11
              (BCType.Interior, BCType.Interior, BCType.Interior, BCType.Interior, BCType.Wall, BCType.PressureFarField)  # b12
              ]

        adj = [((6, 3), (0, 1), 1, True),
               ((0, 2), (0, 0), 0, False),
               ((1, 4), (0, 3), 1, False),
               ((0, 4), (0, 0), 0, False),
               ((0, 0), (0, 5), 1, False),
               ((0, 6), (3, 5), 0, False),

               ((0, 0), (1, 1), 1, False),
               ((1, 2), (0, 0), 0, False),
               ((2, 1), (1, 3), 1, True),
               ((0, 0), (1, 5), 1, False),
               ((1, 6), (4, 5), 0, False),

               ((2, 2), (0, 0), 0, False),
               ((8, 1), (2, 3), 1, True),
               ((2, 4), (0, 0), 0, False),
               ((0, 0), (2, 5), 1, False),
               ((2, 6), (5, 5), 0, False),

               ((9, 3), (3, 1), 1, True),
               ((3, 2), (0, 0), 0, False),
               ((4, 4), (3, 3), 1, False),
               ((3, 4), (0, 0), 0, False),
               ((3, 6), (0, 0), 0, False),

               ((12, 3), (4, 1), 1, True),
               ((4, 2), (0, 0), 0, False),
               ((5, 1), (4, 3), 1, True),
               ((4, 6), (0, 0), 0, False),

               ((5, 2), (0, 0), 0, False),
               ((11, 1), (5, 3), 1, True),
               ((5, 4), (0, 0), 0, False),
               ((5, 6), (0, 0), 0, False),

               ((0, 0), (6, 1), 1, False),
               ((6, 2), (0, 0), 0, False),
               ((6, 4), (7, 2), 0, True),
               ((0, 0), (6, 5), 1, False),
               ((6, 6), (9, 5), 0, False),

               ((8, 2), (7, 1), 1, False),
               ((0, 0), (7, 3), 1, False),
               ((7, 4), (0, 0), 0, False),
               ((0, 0), (7, 5), 1, False),
               ((7, 6), (10, 5), 0, False),

               ((0, 0), (8, 3), 1, False),
               ((8, 4), (0, 0), 0, False),
               ((0, 0), (8, 5), 1, False),
               ((8, 6), (11, 5), 0, False),

               ((12, 2), (9, 1), 1, False),
               ((9, 2), (0, 0), 0, False),
               ((9, 4), (10, 2), 0, True),
               ((9, 6), (0, 0), 0, False),

               ((11, 2), (10, 1), 1, False),
               ((12, 4), (10, 3), 1, False),
               ((10, 4), (0, 0), 0, False),
               ((10, 6), (0, 0), 0, False),

               ((12, 1), (11, 3), 1, True),
               ((11, 4), (0, 0), 0, False),
               ((11, 6), (0, 0), 0, False),

               ((0, 0), (12, 5), 1, False),
               ((12, 6), (0, 0), 0, False)]

        '''构建MSH文件'''
        msh = XF_MSH.from_str3d_multi(blk, bc, adj)
        msh.save(fn)


class HorizontalStablizer(object):
    def __init__(self):
        self.surf = None
        self.tail = None

    @classmethod
    def from_section_param(cls, foil, zoff, cl, swp_bk, tws, dih, ptws, rfy, tkf, pan_dir):
        """
        Construct the vertical stablizer from section geometrical parameters.
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
        Construct the vertical stablizer from section geometrical parameters.
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
