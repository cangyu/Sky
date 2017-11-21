import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from scipy.integrate import romberg
from scipy.optimize import root
from grid import chebshev_dist_multi
from wing import Wing
from aircraft.baseline import WingPlanform


class BWBPlanform(WingPlanform):
    def __init__(self, *args, **kwargs):
        """
        Parametric modeling for aircraft in BWB(Blended Wing Body) configuration.
        :param args: Geometric parameters.
        :param kwargs: Options.
        """

        '''Basic geometric parameters'''
        assert len(args) == 7
        self.Cr = args[0]  # Chord length at root
        self.Cm = args[1]  # Chord length at mid
        self.Ct = args[2]  # Chord length at tip
        self.Bm = args[3]  # Width of the inner part
        self.Bt = args[4]  # Half span
        self.Am = args[5]  # Mean sweep-back of the inner part
        self.At = args[6]  # Mean sweep-back at tip.

        '''Calculate pivots on each curve'''
        front_pnt = np.zeros((3, 3))
        tail_pnt = np.zeros((3, 3))
        tail_pnt[0][0] = self.Cr
        front_pnt[1][0] = self.Bm * math.tan(math.radians(self.Am))
        front_pnt[1][2] = self.Bm
        tail_pnt[1][0] = front_pnt[1][0] + self.Cm
        tail_pnt[1][2] = self.Bm
        front_pnt[2][0] = front_pnt[1][0] + (self.Bt - self.Bm) * math.tan(math.radians(self.At))
        front_pnt[2][2] = self.Bt
        tail_pnt[2][0] = front_pnt[2][0] + self.Ct
        tail_pnt[2][2] = self.Bt

        '''Build interpolated functions'''
        self.u = np.array([0, self.Bm / self.Bt, 1.0])
        self.zw = np.array([0, self.Bm, self.Bt])
        self.xf = make_interp_spline(self.u, front_pnt[:, 0], 3, bc_type=([(1, 0)], [(2, 0)]))
        self.xt = make_interp_spline(self.u, tail_pnt[:, 0], 3, bc_type=([(1, 0)], [(2, 0)]))

    def x_front(self, u):
        return self.xf(u)

    def x_tail(self, u):
        return self.xt(u)

    def z(self, u):
        return u * self.Bt

    def __repr__(self):
        a0 = self.area
        a1 = self.mean_aerodynamic_chord
        a2 = self.span
        a3 = self.root_chord_len
        a4 = self.tip_chord_len
        a5 = a3 / a4
        a6 = a2 ** 2 / a0

        ret = "Blended-Wing-Body Parametric Modeling:\n"
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


def solve_mid_param(*args, **kwargs):
    """
    Calculate the middle chord length and inner width with constant Area and MAC.
    :return: Approximate 'b_mid' and 'c_mid'.
    """

    area2 = args[0]  # Half of the area of the wing.
    mac = args[1]  # Mean aerodynamic chord length of the wing.
    span2 = args[2]  # Half of the width of the span.
    c_root = args[3]  # Length of the root chord.
    c_tip = args[4]  # Length of the tip chord.
    alpha_mid = args[5]  # Averaged sweep back angle of the inner wing.
    alpha_tip = args[6]  # Averaged sweep back angle of the outer wing.

    '''Initial guess of b_mid and c_mid'''
    init_guess = kwargs['init'] if 'init' in kwargs else np.array([0.3 * span2, 0.5 * c_root])

    '''Local constants'''
    tg1 = math.tan(math.radians(alpha_mid))
    tg2 = math.tan(math.radians(alpha_tip))

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

    ans = root(f, init_guess)
    return ans.x


def main():
    '''Wing Planform'''
    c_root = 11
    c_mid = 5.6
    c_tip = c_mid / 2.6
    b_mid = 6.9
    b_tip = 21
    alpha_mid = 32
    alpha_tip = 28

    '''Profile distribution'''
    inner_sec_num = 6
    outer_sec_num = 9
    u_mid = b_mid / b_tip
    u_dist = chebshev_dist_multi([0, u_mid, 1], [inner_sec_num, outer_sec_num])
    sec_num = len(u_dist)

    frm = BWBPlanform(c_root, c_mid, c_tip, b_mid, b_tip, alpha_mid, alpha_tip)
    print(frm)
    frm.show(u_dist)

    '''Profile details'''
    foil = ['SC(2)-0414', 'SC(2)-0414', 'SC(2)-0612', 'SC(2)-0712', 'SC(2)-0710', 'SC(2)-0710', 'SC(2)-0710',
            'SC(2)-1010', 'SC(2)-1010', 'SC(2)-1006', 'SC(2)-0706', 'SC(2)-0706', 'SC(2)-0606', 'SC(2)-0406']
    z_offset = list(map(frm.z, u_dist))
    length = list(map(frm.chord_len, u_dist))
    sweep_back = list(map(lambda u: math.degrees(math.atan2(frm.x_front(u), frm.z(u))), u_dist))
    twist = np.zeros(sec_num)
    dihedral = np.zeros(sec_num)
    twist_pos = np.full(sec_num, 0.25)
    y_ref = np.zeros(sec_num)
    thickness_factor = np.ones(sec_num)

    # '''Record sampling results'''
    # fsp = open('sample_result.txt', 'w')
    # for _spk, param in enumerate(sp):
    #     fsp.write("Case {}:\n".format(_spk))
    #     fsp.write("{:^10} {:^10}\n".format('Profile', 'Twist'))
    #     for l in range(sec_num):
    #         fsp.write("{:^10}{:^10.2f}\n".format(l, param[l]))
    #     fsp.write('\n')
    # fsp.close()

    '''Generate wings'''
    wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)
    model = wg.iges_model()
    model.save('BWB.igs')


if __name__ == '__main__':
    main()
