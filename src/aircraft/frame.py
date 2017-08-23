import math
from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline
from scipy.integrate import romberg
from scipy.optimize import fsolve


class WingFrame(object):
    def __init__(self, xf, xt, yf, yt, z):
        """
        机翼外框描述
        :param xf: 前缘x坐标的参数方程
        :param xt: 后缘x坐标的参数方程
        :param yf: 前缘y坐标的参数方程
        :param yt: 后缘y坐标的参数方程
        :param z: z坐标的参数方程
        """

        self.f_xf = deepcopy(xf)
        self.f_xt = deepcopy(xt)
        self.f_yf = deepcopy(yf)
        self.f_yt = deepcopy(yt)
        self.f_z = deepcopy(z)

    def x_front(self, u):
        return self.f_xf(u)

    def x_tail(self, u):
        return self.f_xt(u)

    def y_front(self, u):
        return self.f_yf(u)

    def y_tail(self, u):
        return self.f_yt(u)

    def z(self, u):
        return self.f_z(u)

    def chord_len(self, u):
        return self.x_tail(u) - self.x_front(u)

    @property
    def area(self):
        """
        Total area of the planar wing. (Not a half)
        :return: The area.
        :rtype: float
        """

        return 2 * math.fabs(romberg(self.chord_len, 0, 1)) * self.z(1)

    @property
    def mean_aerodynamic_chord(self):
        """
        Get the mean aerodynamic chord length of the wing.
        :return: The MAC.
        :rtype: float
        """

        return self.z(1) * math.fabs(romberg(lambda u: self.chord_len(u) ** 2, 0, 1)) / (0.5 * self.area)

    @property
    def span(self):
        return 2 * (self.z(1) - self.z(0))

    @property
    def root_chord_len(self):
        return self.chord_len(0)

    @property
    def tip_chord_len(self):
        return self.chord_len(1)

    def __str__(self):
        ret = "Wing Planar Frame Info:\n"
        ret += "Area: {:.4f} m^2\n".format(self.area)
        ret += "MAC: {:.3f} m\n".format(self.mean_aerodynamic_chord)
        ret += "Span: {:.3f} m\n".format(self.span)
        ret += "Root: {:.3f} m\n".format(self.root_chord_len)
        ret += "Tip: {:.3f} m".format(self.tip_chord_len)
        return ret

    def show(self, n=1000):
        u_dist = np.linspace(0, 1.0, n)
        z = np.empty(n, float)
        xf = np.empty(n, float)
        xt = np.empty(n, float)
        for k in range(n):
            z[k] = self.z(u_dist[k])
            xf[k] = self.x_front(u_dist[k])
            xt[k] = self.x_tail(u_dist[k])

        plt.figure()
        plt.plot(z, xf, label='Front')
        plt.plot(z, xt, label='Tail')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.show()


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

        self.Cr = c_root
        self.Cm = c_mid
        self.Ct = c_tip
        self.Bm = b_mid
        self.Bt = b_tip
        self.Am = alpha_mid
        self.At = alpha_tip

        '''Calculate pivots on each curve'''
        front_pnt = np.empty((3, 3), float)
        tail_pnt = np.empty((3, 3), float)
        front_pnt[0] = np.zeros(3)
        tail_pnt[0] = np.array([self.Cr, 0, 0])
        front_pnt[1] = np.array([self.Bm * math.tan(math.radians(self.Am)), 0, self.Bm])
        tail_pnt[1] = np.array([front_pnt[1][0] + self.Cm, 0, self.Bm])
        front_pnt[2] = np.array([front_pnt[1][0] + (self.Bt - self.Bm) * math.tan(math.radians(self.At)), 0, self.Bt])
        tail_pnt[2] = np.array([front_pnt[2][0] + self.Ct, 0, self.Bt])

        '''Build interpolated functions'''
        u = np.array([0, self.Bm / self.Bt, 1.0])
        z = np.array([0, self.Bm, self.Bt])
        xf = np.array([front_pnt[0][0], front_pnt[1][0], front_pnt[2][0]])
        yf = np.array([front_pnt[0][1], front_pnt[1][1], front_pnt[2][1]])
        xt = np.array([tail_pnt[0][0], tail_pnt[1][0], tail_pnt[2][0]])
        yt = np.array([tail_pnt[0][1], tail_pnt[1][1], tail_pnt[2][1]])

        super(BWBFrame, self).__init__(make_interp_spline(u, xf, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       make_interp_spline(u, xt, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       make_interp_spline(u, yf, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       make_interp_spline(u, yt, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       lambda t: z[1] * t / u[1] if t <= u[1] else z[1] + (z[2] - z[1]) * (t - u[1]) / (u[2] - u[1]))

    @classmethod
    def from_area_mac(cls, area, c_mid, mac, b_mid, b_tip, alpha_mid, alpha_tip):
        """
        Construct the planar frame with constant Area and MAC.
        :param area: Area of the wing. (Not a half)
        :type area: float
        :param c_mid: Length of the middle chord.
        :type c_mid: float
        :param mac: Mean aerodynamic chord length of the wing.
        :type mac: float
        :param b_mid: Width of the inner wing.
        :type b_mid: float
        :param b_tip: Width of the half of the span.
        :type b_tip: float
        :param alpha_mid: Averaged sweep back angle of the inner wing.
        :type alpha_mid: float
        :param alpha_tip: Averaged sweep back angle of the outer wing.
        :type alpha_tip: float
        :return: Constrained frame.
        :rtype: BWBFrame
        """

        '''Calculate pivots on each curve'''
        area2 = area / 2
        p = np.zeros((6, 3))

        p[1][0] = b_mid * math.tan(alpha_mid)
        p[1][2] = b_mid

        p[2][0] = p[1][0] + (b_tip - b_mid) * math.tan(alpha_tip)
        p[5][2] = p[2][2] = b_tip

        p[4] = p[1]
        p[4][0] += c_mid

        u = np.array([0, b_mid / b_tip, 1.0])
        xf = make_interp_spline(u, p[:3][0], 3, bc_type=([(1, 0)], [(2, 0)]))

        def sp(_cr, _ct):
            tc = make_interp_spline(u, [_cr, c_mid, _ct], 3, bc_type=([(1, 0)], [(2, 0)]))
            ts = romberg(lambda _u: tc(_u) - xf(_u), 0, 1)
            tmac = romberg(lambda _u: (tc(_u) - xf(_u)) ** 2, 0, 1) / ts
            return ts - area2, tmac - mac

        taper_ratio = 5
        pinit = np.array([taper_ratio, 1]) * (area / b_tip) / (1 + taper_ratio)
        ans = fsolve(sp, pinit)
        p[3][0] = ans[0]
        p[5][0] = ans[1]

        return BWBFrame(p[3][0], c_mid, p[5][0], b_mid, b_tip, alpha_mid, alpha_tip)
