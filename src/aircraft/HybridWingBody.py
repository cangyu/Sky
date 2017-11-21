import math
import numpy as np
from matplotlib import pylab as plt
from grid import chebshev_dist_multi
from aircraft import WingFrame
from wing import Wing
from iges import Model


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



'''Wing Planform'''
spn = 21
cr = 16.667
ct = 1.38
outer_taper_ratio = 2.5
fl = np.array([2.1, 4, 2.9, 0])
alpha = np.radians([32, 56, 37.5, 28])
tl = np.array([1.5, 3.2, 0, 0])
beta = np.radians([-26, -40, 0, 0])
frm = HWBFrame(spn, cr, ct, fl, alpha, tl, beta, outer_taper=outer_taper_ratio)
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

foil = ['NACA14122',
        'NACA14022',
        'NACA63(4)-221',
        'NACA63(3)-218',
        'NLF(1)-0416',
        'NLF(1)-0414F',
        'SC(2)-0612',
        'SC(2)-0610',
        'SC(2)-0710',
        'SC(2)-0710',
        'SC(2)-0610',
        'SC(2)-0410',
        'NACA64A210',
        'NACA64A010']
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
model.add(wg.surf().to_iges())
model.save('HWB_Wing.igs')

wg.gen_grid('HWB_Wing.msh')
