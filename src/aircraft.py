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
