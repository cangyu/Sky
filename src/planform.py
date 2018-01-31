#!/usr/bin/env python
# -*- coding: utf-8 -*-

from abc import abstractmethod, ABCMeta
import math
import numpy as np
from scipy.integrate import romberg
from scipy.misc import derivative
from spacing import uniform
from misc import pnt_pan, share, normalize
from nurbs import LocalCubicInterpolatedCrv, point_inverse, Spline
from settings import z_axis_positive


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
    def y_tail(self, u):
        """
        Calculate the Y-coordinate on trailing edge.
        :param u: Relative position parameter.
        :type u: float
        :return: Y-coordinate on trailing edge.
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

    def x_025(self, u):
        return share(0.25, self.x_front(u), self.x_tail(u))

    def swp_025(self, u, delta=1e-3):
        tangent = derivative(self.x_025, u, dx=delta)
        return math.degrees(math.atan2(tangent, self.span / 2))

    def chord_len(self, u):
        return self.x_tail(u) - self.x_front(u)

    def area_between(self, u_begin, u_end):
        spn2 = self.span / 2
        return romberg(self.chord_len, u_begin, u_end) * spn2

    @property
    def area(self):
        """
        Total area of the planar wing. (Not a half)
        :return: The area.
        :rtype: float
        """

        return 2 * self.area_between(0, 1)

    @property
    def mean_aerodynamic_chord(self):
        """
        Get the mean aerodynamic chord length of the wing.
        :return: The MAC.
        :rtype: float
        """

        spn2 = self.span / 2
        ret = romberg(lambda u: self.chord_len(u) ** 2, 0, 1)
        ret *= spn2
        ret /= 0.5 * self.area
        return ret

    @property
    def span(self):
        return 2 * (self.z(1) - self.z(0))

    @property
    def half_span(self):
        return 0.5 * self.span

    @property
    def root_chord_len(self):
        return self.chord_len(0)

    @property
    def tip_chord_len(self):
        return self.chord_len(1)

    def __repr__(self):
        return 'Wing Planform'

    def pic(self, *args, **kwargs):
        ax = args[0]
        n = args[1] if len(args) == 2 else 100

        u = uniform(n)
        spn2 = self.span / 2

        z = np.array([self.z(t) for t in u])
        leading_x = np.array([self.x_front(t) for t in u])
        trailing_x = np.array([self.x_tail(t) for t in u])

        if 'direction' in kwargs and kwargs['direction'] == 'vertical':
            ax.plot(leading_x, z, label='Leading')
            ax.plot(trailing_x, z, label='Trailing')
            if 'u' in kwargs:
                u_pos = np.copy(kwargs['u'])
                for lu in u_pos:
                    zp = lu * spn2
                    ax.plot([self.x_front(lu), self.x_tail(lu)], [zp, zp], '--')
        elif 'direction' not in kwargs or kwargs['direction'] == 'horizontal':
            ax.plot(z, leading_x, label='Leading')
            ax.plot(z, trailing_x, label='Trailing')
            if 'u' in kwargs:
                u_pos = np.copy(kwargs['u'])
                for k, lu in enumerate(u_pos):
                    zp = lu * spn2
                    xf = self.x_front(lu)
                    xt = self.x_tail(lu)
                    ax.plot([zp, zp], [xf, xt], '--')
                    ax.text(zp, share(0.5, xf, xt), str(k))

            ax.invert_yaxis()
        else:
            raise AttributeError('invalid direction')

        ax.set_aspect('equal')
        ax.legend()


class HWBInnerStraightPlanform(WingPlanform):
    def __init__(self, *args, **kwargs):
        """
        Wing Planform for aircraft under Hybrid-Wing-Body configuration.
        :param args: Geometric parameters describing the shape.
                        wing_root_len, wing_tip_len, wing_spn2,
                        wing_leading_inner_delta, wing_leading_middle_delta, wing_leading_outer_sweep,
                        wing_trailing_inner_delta, wing_trailing_outer_spn, wing_trailing_outer_sweep
        :param kwargs: Options.
        """

        cr, ct, spn2 = args[0:3]
        leading_inner_delta, leading_middle_delta, leading_outer_swp = args[3:6]
        trailing_inner_delta, trailing_outer_spn, trailing_outer_swp = args[6:9]

        leading_seg_length = np.empty(3)
        leading_seg_length[0] = leading_inner_delta[0]
        leading_seg_length[1] = leading_middle_delta[0]
        leading_seg_length[2] = spn2 - (leading_seg_length[0] + leading_seg_length[1])

        trailing_seg_length = np.empty(3)
        trailing_seg_length[0] = trailing_inner_delta[0]
        trailing_seg_length[2] = trailing_outer_spn
        trailing_seg_length[1] = spn2 - (trailing_seg_length[0] + trailing_seg_length[2])

        leading_seg_theta = np.empty(3)
        leading_seg_theta[0] = math.atan2(leading_inner_delta[1], leading_inner_delta[0])
        leading_seg_theta[1] = math.atan2(leading_middle_delta[1], leading_middle_delta[0])
        leading_seg_theta[2] = math.radians(leading_outer_swp)

        trailing_seg_theta = np.empty(3)
        trailing_seg_theta[0] = math.atan2(trailing_inner_delta[1], trailing_inner_delta[0])
        trailing_seg_theta[2] = math.radians(trailing_outer_swp)
        leading_dx = sum([leading_seg_length[i] * math.tan(leading_seg_theta[i]) for i in range(3)])
        dx = leading_dx + ct - math.tan(trailing_seg_theta[2]) * trailing_seg_length[2] - (cr + trailing_inner_delta[1])
        dz = trailing_seg_length[1]
        trailing_seg_theta[1] = math.atan2(dx, dz)

        desc_shape = (4, 3)
        leading_pnt = np.empty(desc_shape, float)
        leading_tangent = np.empty(desc_shape, float)
        trailing_pnt = np.empty(desc_shape, float)
        trailing_tangent = np.empty(desc_shape, float)

        leading_pnt[0] = (0, 0, 0) if 'origin' not in kwargs else kwargs['origin']
        for i in range(1, 4):
            delta = (math.tan(leading_seg_theta[i - 1]) * leading_seg_length[i - 1], 0, leading_seg_length[i - 1])
            leading_pnt[i] = pnt_pan(leading_pnt[i - 1], delta)

        trailing_pnt[0] = pnt_pan(leading_pnt[0], (cr, 0, 0))
        for i in range(1, 4):
            delta = (math.tan(trailing_seg_theta[i - 1]) * trailing_seg_length[i - 1], 0, trailing_seg_length[i - 1])
            trailing_pnt[i] = pnt_pan(trailing_pnt[i - 1], delta)

        leading_tangent[0] = leading_tangent[1] = (math.sin(leading_seg_theta[0]), 0, math.cos(leading_seg_theta[0]))
        leading_tangent[2] = leading_tangent[3] = (math.sin(leading_seg_theta[2]), 0, math.cos(leading_seg_theta[2]))

        trailing_tangent[0] = trailing_tangent[1] = (math.sin(trailing_seg_theta[0]), 0, math.cos(trailing_seg_theta[0]))
        trailing_tangent[2] = trailing_tangent[3] = (math.sin(trailing_seg_theta[2]), 0, math.cos(trailing_seg_theta[2]))

        self.Span2 = spn2
        self.LeadingCrv = LocalCubicInterpolatedCrv(leading_pnt, leading_tangent)
        self.TrailingCrv = LocalCubicInterpolatedCrv(trailing_pnt, trailing_tangent)

    def z(self, u):
        return u * self.Span2

    def x_front(self, u):
        zp = self.z(u)
        lu = point_inverse(self.LeadingCrv, zp, 2)
        return self.LeadingCrv(lu)[0]

    def x_tail(self, u):
        zp = self.z(u)
        lu = point_inverse(self.TrailingCrv, zp, 2)
        return self.TrailingCrv(lu)[0]

    def y_front(self, u):
        return 0

    def y_tail(self, u):
        return 0


class HWBNoseBluntPlanform(WingPlanform):
    def __init__(self, *args, **kwargs):
        root_chord, mid_chord, tip_chord, p1, p2, p3, p4, p5, p6, p7, p8, p9 = args
        tension = np.ones(9) if 'tension' not in kwargs else kwargs['tension']
        assert len(tension) == 9

        self.p = np.zeros((13, 3))
        self.t = np.zeros((9, 3))

        if 'ref' in kwargs:
            tmp = kwargs['ref']
            for i in range(3):
                self.p[0][i] = tmp[i]
        for i in range(1, 13):
            self.p[i][1] = self.p[0][1]

        def assign_helper(dst, src):
            dst[0] = src[0]
            dst[2] = src[1]

        assign_helper(self.p[1], p1)
        assign_helper(self.p[2], p2)
        assign_helper(self.p[3], p3)
        assign_helper(self.p[4], p4)
        assign_helper(self.p[5], p5)
        assign_helper(self.p[6], p6)
        assign_helper(self.p[7], p7)
        self.p[8] = pnt_pan(self.p[0], (root_chord, 0, 0))
        assign_helper(self.p[9], p8)
        self.p[10] = pnt_pan(self.p[6], (mid_chord, 0, 0))
        assign_helper(self.p[11], p9)
        self.p[12] = pnt_pan(self.p[7], (tip_chord, 0, 0))

        self.t[0] = z_axis_positive
        self.t[2] = self.t[1] = normalize(self.p[3] - self.p[2])
        self.t[4] = self.t[3] = normalize(self.p[7] - self.p[6])
        self.t[6] = self.t[5] = normalize(self.p[9] - self.p[8])
        self.t[8] = self.t[7] = normalize(self.p[12] - self.p[11])

        for i in range(9):
            self.t[i] *= tension[i]

        self.seg1 = Spline(self.p[:3], p=3, bc=([(1, self.t[0])], [(1, self.t[1])]))
        self.seg2 = Spline(self.p[3:7], p=3, bc=([(1, self.t[2])], [(1, self.t[3])]))
        self.seg3 = Spline(self.p[9:12], p=3, bc=([(1, self.t[6])], [(1, self.t[7])]))

    @property
    def span(self):
        return 2 * (self.p[7][2] - self.p[0][2])

    def z(self, u):
        return u * self.half_span

    def x_front(self, u):
        z = self.z(u)

        if z < self.p[0][2]:
            raise ValueError('out of boundary')
        elif z < self.p[2][2]:
            lu = point_inverse(self.seg1, z, 2)
            return self.seg1(lu)[0]
        elif z < self.p[3][2]:
            ratio = (z - self.p[2][2]) / (self.p[3][2] - self.p[2][2])
            return share(ratio, self.p[2][0], self.p[3][0])
        elif z < self.p[6][2]:
            lu = point_inverse(self.seg2, z, 2)
            return self.seg2(lu)[0]
        elif z <= self.p[7][2]:
            ratio = (z - self.p[6][2]) / (self.p[7][2] - self.p[6][2])
            return share(ratio, self.p[6][0], self.p[7][0])
        else:
            raise ValueError('out of boundary')

    def x_tail(self, u):
        z = self.z(u)

        if z < self.p[8][2]:
            raise ValueError('out of boundary')
        elif z < self.p[9][2]:
            ratio = (z - self.p[8][2]) / (self.p[9][2] - self.p[8][2])
            return share(ratio, self.p[8][0], self.p[9][0])
        elif z < self.p[11][2]:
            lu = point_inverse(self.seg3, z, 3)
            return self.seg3(lu)[0]
        elif z <= self.p[12][2]:
            ratio = (z - self.p[11][2]) / (self.p[12][2] - self.p[11][2])
            return share(ratio, self.p[11][0], self.p[12][0])
        else:
            raise ValueError('out of boundary')

    def y_front(self, u):
        return self.p[0][1]

    def y_tail(self, u):
        return self.p[0][1]
