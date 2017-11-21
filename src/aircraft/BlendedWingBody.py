import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from scipy.integrate import romberg
from scipy.optimize import root
from grid import chebshev_dist_multi
from wing import Wing
from iges import Model
from nurbs import point_inverse, LocalCubicInterpolatedCrv, Line, Coons
from aircraft.Baseline import WingPlanform


class BWBPlanform1(WingPlanform):
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


class BWBPlanform2(WingPlanform):
    def __init__(self, spn, cr, ct, fl, alpha, tl, beta, outer_taper=2.5):
        """
        Parametric wing Planform for BWB configuration.
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
            fp[i + 1][0] += fl[i] * math.tan(alpha[i])
            fp[i + 1][2] += fl[i]
        ftv = np.array([[0, 0, 1],
                        [math.sin(alpha[1]), 0, math.cos(alpha[1])],
                        [math.sin(alpha[1]), 0, math.cos(alpha[1])],
                        [math.sin(alpha[3]), 0, math.cos(alpha[3])],
                        [math.sin(alpha[3]), 0, math.cos(alpha[3])]])

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
            tp[i + 1][0] += tl[i] * math.tan(beta[i])
            tp[i + 1][2] += tl[i]

        tp[-1] = fp[-1]
        tp[-1][0] += ct

        tp[-2] = fp[-2]
        tp[-2][0] += outer_taper * ct

        beta[-1] = math.atan2(tp[-1][0] - tp[-2][0], tl[-1])
        beta[-2] = math.atan2(tp[-2][0] - tp[-3][0], tl[-2])

        ttv = np.array([[0, 0, 1],
                        [math.sin(beta[1]), 0, math.cos(beta[1])],
                        [math.sin(beta[1]), 0, math.cos(beta[1])],
                        [math.sin(beta[3]), 0, math.cos(beta[3])],
                        [math.sin(beta[3]), 0, math.cos(beta[3])]])

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

    def __repr__(self):
        a0 = self.area
        a1 = self.mean_aerodynamic_chord
        a2 = self.span
        a3 = self.root_chord_len
        a4 = self.tip_chord_len
        a5 = a3 / a4
        a6 = a2 ** 2 / a0

        ret = "BWB2 Wing Planform Info:\n"
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


def planform1():
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

    frm = BWBPlanform1(c_root, c_mid, c_tip, b_mid, b_tip, alpha_mid, alpha_tip)
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


def planform2():
    spn = 21
    cr = 16.667
    ct = 1.38
    outer_taper_ratio = 2.5
    fl = np.array([2.1, 4, 2.9, 0])
    alpha = np.radians([32, 56, 37.5, 28])
    tl = np.array([1.5, 3.2, 0, 0])
    beta = np.radians([-26, -40, 0, 0])
    frm = BWBPlanform2(spn, cr, ct, fl, alpha, tl, beta, outer_taper=outer_taper_ratio)
    print(frm)

    '''Profile distribution'''
    u = [0]
    tmp = 0
    for l in fl:
        tmp += l
        u.append(tmp)
    tmp = 0
    for l in tl:
        tmp += l
        u.append(tmp)
    u = np.unique(u) / spn
    u_dist = chebshev_dist_multi([0, u[1], u[2], u[3], u[4], u[5], u[6]], [3, 2, 3, 3, 3, 5])
    frm.show(u_dist)

    '''Profile details'''
    sec_num = len(u_dist)
    front_swp = np.zeros(sec_num)
    for i in range(sec_num):
        tg = frm.front_crv(u_dist[i], 1)
        front_swp[i] = math.degrees(math.atan2(tg[0], tg[2]))

    tc_3d = np.array([22, 22, 21, 19, 16, 14, 12, 11, 10, 8, 8, 8, 8, 8], float)
    cl_3d = np.array([0.08, 0.10, 0.14, 0.18, 0.27, 0.38, 0.42, 0.45, 0.48, 0.47, 0.44, 0.29, 0.1, 0])
    cl_2d = 1.1 * np.copy(cl_3d) / math.cos(math.radians(28)) ** 2
    print(cl_2d)

    foil = ['NACA14122', 'NACA14022', 'NACA63(4)-221', 'NACA63(3)-218', 'NLF(1)-0416', 'NLF(1)-0414F', 'SC(2)-0612',
            'SC(2)-0610', 'SC(2)-0710', 'SC(2)-0710', 'SC(2)-0610', 'SC(2)-0410', 'NACA64A210', 'NACA64A010']
    z_offset = u_dist * spn
    length = list(map(lambda _u: frm.chord_len(_u), u_dist))
    sweep_back = list(map(lambda _u: math.degrees(math.atan2(frm.x_front(_u), frm.z(_u))), u_dist))
    twist = np.zeros(sec_num)
    dihedral = np.zeros(sec_num)
    twist_pos = np.ones(sec_num)
    y_ref = np.zeros(sec_num)
    thickness_factor = np.ones(sec_num)

    '''Show distribution'''
    f, ax_arr = plt.subplots(2, sharex=True)
    ax_arr[0].plot(u_dist * spn, tc_3d / 100)
    ax_arr[0].set_ylim([0, tc_3d[0] / 100 * 1.1])
    ax_arr[0].set_title('t/c in span-wise')
    ax_arr[1].plot(u_dist * spn, cl_3d)
    ax_arr[1].set_title('Cl in span-wise')
    plt.show()

    '''Initial grid'''
    wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)
    model = Model()
    sf = wg.surf
    model.add(sf.to_iges())
    model.save('HWB_Wing.igs')
    # wg.gen_grid('HWB_Wing.msh')


if __name__ == '__main__':
    # planform1()
    planform2()
