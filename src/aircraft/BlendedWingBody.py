import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.interpolate import make_interp_spline
from scipy.integrate import romberg
from scipy.optimize import root
from grid import chebshev_dist_multi, uniform
from wing import Wing
from iges import Model, Entity110, Entity116
from nurbs import point_inverse, LocalCubicInterpolatedCrv, Line, Coons
from aircraft.Baseline import WingPlanform
from misc import share


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


def tangent_on_crv(pos, crv):
    tg = crv(pos, 1)
    return math.degrees(math.atan2(tg[0], tg[2]))


def planform2():
    """
    Build a BWB with more detailed control.
    :return: None
    """

    '''Planform parameters'''
    spn = 21
    cr = 28
    ct = 1.5
    fl = np.array([1.2, 2.9, 3.9, 0])  # 最后一个由程序自动计算
    alpha = np.radians([60.0, 71.0, 59.8, 26])
    tl = np.array([4.9, 0.92, 0, 0])  # 后两个由程序自动计算
    beta = np.radians([-23.8, -65.0, 0, 0])  # 后两个由程序自动计算
    frm = BWBPlanform2(spn, cr, ct, fl, alpha, tl, beta, outer_taper=4.0)
    # print(frm)

    w = 120 * 1e3 * 9.8
    a = 299.5
    rho = 0.4135
    s = 360.31
    for Ma in [0.65, 0.7, 0.75, 0.8, 0.85]:
        p_inf = 0.5 * rho * (Ma * a) ** 2
        print('Design Cl when Ma = {:>4.2f}: {:>4.2f}'.format(Ma, w / (p_inf * s)))

    '''Profile distribution'''
    z_pos = np.array([0, 1.5, 3, 5, 6.5, 8, 10, 13, 17, 20, 21])
    u_pos = np.array([z / spn for z in z_pos])
    n = len(u_pos)

    '''Profile details'''
    front_sweep = np.array([tangent_on_crv(u, frm.front_crv) for u in u_pos])
    chord_len = np.array([frm.chord_len(u) for u in u_pos])
    tc = np.array([16.0, 16.0, 16.0, 14.0, 12.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0])
    height = chord_len * tc / 100
    cl = np.array([0.100, 0.15, 0.21, 0.30, 0.35, 0.40, 0.39, 0.37, 0.28, 0.13, 0.0])  # for airfoil
    cl_airfoil = np.array([1.1 * cl[i] / math.cos(math.radians(front_sweep[i])) ** 2 for i in range(n)])
    foil = ['NACA05116', 'NACA15116', 'NACA15116', 'NACA15114', 'BL0', 'BL1256',
            'SC(2)-0410', 'SC(2)-0410', 'SC(2)-0410', 'SC(2)-0410', 'SC(2)-0010']
    twist = np.array([-0.028, -0.025, 0.271, 0.790, 1.155, 0.495,
                      0.160, 0.068, -0.362, -1.091, 0])

    sweep_back = np.array([math.degrees(math.atan2(frm.x_front(u), frm.z(u))) for u in u_pos])
    dihedral = np.array([math.degrees(math.atan2((height[0] - height[i]) / 2, z_pos[i])) for i in range(n)])
    twist_pos = np.ones(n)
    y_ref = np.zeros(n)
    thickness_factor = np.ones(n)

    print('\n{:^4}{:>8}{:>8}{:>14}{:>14}{:>12}{:>16}{:>12}'.format('idx', 'z/m', 'pos/%', 'chord len/m', 'max height/m', 'profile cl', 'front sweep/deg', 'airfoil cl'))
    for i in range(n):
        tmp = '{:^4}'.format(i)
        tmp += '{:>8.3f}{:>8.2f}'.format(z_pos[i], 100 * u_pos[i])
        tmp += '{:>10.3f}    {:>10.3f}    '.format(chord_len[i], height[i])
        tmp += '{:>10.3f}  {:>12.2f}    {:>10.3f}  '.format(cl[i], front_sweep[i], cl_airfoil[i])
        print(tmp)

    '''CAD Geom'''
    model = Model()
    wg = Wing.from_geom_desc(foil, chord_len, thickness_factor, z_pos, sweep_back, twist, twist_pos, dihedral, y_ref)
    for k, elem in enumerate(wg.profile):
        crv = elem.crv
        if k == 4:
            crv.pan([0, -0.02, 0])
        if k >= 5:
            crv.pan([0, 0.08 + (z_pos[k] - z_pos[4]) * math.tan(math.radians(1)), 0])
        model.add(crv.to_iges())

    # ll = 1200
    # hh = 300
    # ww = 1000
    # far_pnt = np.array([[-ll, hh, 0], [ll, hh, 0], [ll, -hh, 0], [-ll, -hh, 0],
    #                     [-ll, hh, ww], [ll, hh, ww], [ll, -hh, ww], [-ll, -hh, ww]])
    # far_wire = np.array([[0, 4], [4, 7], [7, 3], [3, 0],
    #                      [0, 1], [4, 5], [7, 6], [3, 2],
    #                      [1, 5], [5, 6], [6, 2], [2, 1]])
    #
    # for pnt in far_pnt:
    #     model.add(Entity116(pnt))
    # for line in far_wire:
    #     model.add(Entity110(far_pnt[line[0]], far_pnt[line[1]]))

    model.save('BWB_Wing.igs')

    '''Plot profile info'''
    fig = plt.figure()
    gs = gridspec.GridSpec(4, 2)
    planform_plt = fig.add_subplot(gs[:, 0])
    tc_plt = fig.add_subplot(gs[0, 1])
    cl_plt = fig.add_subplot(gs[1, 1])
    twist_plt = fig.add_subplot(gs[2, 1])
    len_height_plt = fig.add_subplot(gs[3, 1])

    # leading and trailing edges
    fp = np.array([frm.front_crv(p) for p in uniform(1000)])
    tp = np.array([frm.tail_crv(p) for p in uniform(1000)])
    planform_plt.plot(fp[:, 2], fp[:, 0], label='Leading Edge')
    planform_plt.plot(tp[:, 2], tp[:, 0], label='Trailing Edge')
    planform_plt.invert_yaxis()
    planform_plt.set_aspect('equal')

    # add chord on each profile
    fu = np.array([point_inverse(frm.front_crv, pos, 2) for pos in z_pos])
    tu = np.array([point_inverse(frm.tail_crv, pos, 2) for pos in z_pos])
    for k in range(n):
        tfx = frm.front_crv(fu[k])[0]
        ttx = frm.tail_crv(tu[k])[0]
        z = z_pos[k]
        planform_plt.plot([z, z], [tfx, ttx], '--', label='Chord {} {} at {:.2f}%'.format(k, foil[k], 100 * u_pos[k]))
        planform_plt.text(z, (tfx + ttx) / 2, str(k))

    planform_plt.legend()
    planform_plt.set_title('BWB Wing Planform')

    # t/c and height distribution
    tc_plt.plot(u_pos, tc, label='t/c')
    tc_plt.plot(u_pos, height, label='camber height')
    tc_plt.set_title('Thickness and height distribution')
    tc_plt.grid(True)
    tc_plt.legend()
    for i in range(n):
        u = u_pos[i]
        tc_plt.text(u, tc[i], '{:.1f}'.format(tc[i]))
        tc_plt.plot([u, u], [0, tc[i]], '--')
        tc_plt.text(u, height[i], '{:.1f}'.format(height[i]))

    # cl distribution
    cl_plt.plot(u_pos, cl, label='Wing Profile')
    cl_plt.plot(u_pos, cl_airfoil, label='Airfoil Desired')
    cl_plt.legend()
    cl_plt.set_title('Lift coefficient distribution')
    cl_plt.grid(True)
    for i in range(n):
        u = u_pos[i]
        cl_plt.text(u, cl[i], '{:.2f}'.format(cl[i]))
        cl_plt.text(u, cl_airfoil[i], '{:.2f}'.format(cl_airfoil[i]))

    # twist distribution
    twist_plt.plot(u_pos, twist, '-^')
    twist_plt.set_title('Twist distribution')
    twist_plt.grid(True)
    for i in range(n):
        u = u_pos[i]
        twist_plt.text(u, twist[i], '{:.2f}'.format(twist[i]))

    # chord length and height
    len_height_plt.bar([i + 1 for i in range(n)], chord_len)
    len_height_plt.set_title('Chord length')
    len_height_plt.grid(True)
    for i in range(n):
        len_height_plt.text(i + 1, chord_len[i], '{:.1f}'.format(chord_len[i]))

    fig.tight_layout()
    plt.show()


def simple_wing():
    span = 42
    # area = 220
    # mac = 6.52
    front_sweep = 28
    cr = 9.73
    ct = 0.75

    foil = ['NACA0012', 'NACA0010', 'NACA0010', 'NACA0008', 'NACA0006']
    u = np.array([0, 0.25, 0.5, 0.75, 1])
    n = len(u)

    z_offset = u * span
    length = np.array([share(_u, cr, ct) for _u in u])
    sweep_back = np.array([front_sweep] * n)
    twist = np.zeros(n)
    dihedral = np.zeros(n)
    twist_pos = np.ones(n)
    y_ref = np.zeros(n)
    thickness_factor = np.ones(n)

    wg = Wing.from_geom_desc(foil, length, thickness_factor, z_offset, sweep_back, twist, twist_pos, dihedral, y_ref)
    model = Model()
    model.add(wg.surf.to_iges())
    model.add(wg.tailing_up.to_iges())
    model.add(wg.tailing_down.to_iges())
    model.save('StraightWing_{}_{}_{}_{}.igs'.format(span, front_sweep, cr, ct))


if __name__ == '__main__':
    # planform1()
    planform2()
    # simple_wing()
