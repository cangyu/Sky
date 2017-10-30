from abc import ABCMeta, abstractmethod
from math import sin, cos, tan, radians, fabs, atan2

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import romberg
from scipy.interpolate import make_interp_spline
from scipy.optimize import root

from nurbs import Coons
from nurbs import LocalCubicInterpolatedCrv, Line, point_inverse


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
