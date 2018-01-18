from abc import abstractmethod, ABCMeta
from copy import deepcopy
import math
import numpy as np
from numpy.linalg import norm
from scipy.integrate import romberg
from scipy.misc import derivative
from grid import LinearTFI2D, LinearTFI3D, Laplace2D, ThomasMiddlecoff2D
from grid import Plot3D, Plot3DBlock, uniform
from grid import hyperbolic_tangent, single_exponential, double_exponential
from iges import Model, Entity116, Entity110
from misc import pnt_dist, read_airfoil_pts, pnt_pan, share, normalize
from nurbs import Crv, Line, Spline, ConicArc
from nurbs import GlobalInterpolatedCrv, LocalCubicInterpolatedCrv
from nurbs import Surf, Skinned, RuledSurf, point_inverse
from settings import AIRFOIL_LIST
from rotation import pnt_rotate

global_origin = np.array([0., 0., 0.])
x_axis_positive = np.array([1., 0, 0])
x_axis_negative = np.array([-1., 0, 0])
y_axis_positive = np.array([0, 1., 0])
y_axis_negative = np.array([0, -1., 0])
z_axis_positive = np.array([0, 0, 1.])
z_axis_negative = np.array([0, 0, -1.])


class WingProfile(Airfoil):
    def __init__(self, foil, ends):
        """
        3D profile at certain position.
        :param foil: 2D Airfoil.
        :type foil: Airfoil
        :param ends: Starting and ending points of the profile.
        """

        assert type(foil) is Airfoil  # avoid derived objects
        original_chord_len = foil.chord_len
        assert not math.isclose(original_chord_len, 0)
        super(WingProfile, self).__init__(foil)

        new_chord_len = pnt_dist(ends[0], ends[1])
        assert not math.isclose(new_chord_len, 0)
        assert math.isclose(ends[0][2], ends[1][2])
        self.ending = np.copy(ends)

        '''Stretch and Z offset'''
        scale_ratio = new_chord_len / original_chord_len
        z_off = self.ending[0][2]
        for i in range(self.pnt_num):
            self.pts[i][0] *= scale_ratio
            self.pts[i][1] *= scale_ratio
            self.pts[i][2] = z_off

        '''Rotate around trailing'''
        ref = share(0.5, self.pts[0], self.pts[-1])
        delta = self.ending[0] - self.ending[1]
        ang = math.degrees(math.atan2(delta[1], delta[0])) - 180
        self.pts = pnt_rotate(ref, z_axis_positive, ang, self.pts)

        '''Move to ends[1]'''
        delta = self.ending[-1] - share(0.5, self.pts[0], self.pts[-1])
        for i in range(self.pnt_num):
            self.pts[i] += delta

    @property
    def front(self):
        return self.ending[0]

    @property
    def tail(self):
        return self.ending[-1]

    @classmethod
    def from_profile_param(cls, wpp):
        """
        Construct the profile from parameters in higher level.
        :param wpp: Wing profile description object.
        :type wpp: WingProfileParam
        :return: Target profile.
        :rtype: WingProfile
        """

        theta = math.radians(wpp.twist_ang)
        pan_dir1 = np.array([-math.cos(theta), math.sin(theta), 0])
        pan_dir2 = -pan_dir1
        len1 = wpp.chord_len * wpp.twist_ref
        len2 = wpp.chord_len - len1

        ending = np.empty((2, 3), float)
        ending[0] = pnt_pan(wpp.twist_center, pan_dir1 * len1)
        ending[1] = pnt_pan(wpp.twist_center, pan_dir2 * len2)

        return cls(wpp.airfoil, ending)


class WingProfileParam(object):
    def __init__(self, *args, **kwargs):
        """
        Intrinsic descriptions for profile in span-wise direction.
        :param args: Elementary parameters.
        :param kwargs: Optional parameters.
        """

        # Airfoil object
        # Chord length, input in meters
        # Angle of twist, input in degrees
        # Spatial position of the center for twisting, also indicates the position of the profile
        # Twist position along the chord

        if len(args) == 1:
            self.airfoil = args[0]
            self.chord_len = 1.0
            self.twist_ang = 0.0
            self.twist_center = np.array([self.chord_len, 0, 0])
            self.twist_ref = 1.0
        elif len(args) == 2:
            self.airfoil = args[0]
            self.chord_len = args[1]
            self.twist_ang = 0.0
            self.twist_center = np.array([self.chord_len, 0, 0])
            self.twist_ref = 1.0
        elif len(args) == 3:
            self.airfoil = args[0]
            self.chord_len = args[1]
            self.twist_ang = args[2]
            self.twist_center = np.array([self.chord_len, 0, 0])
            self.twist_ref = 1.0
        elif len(args) == 4:
            raise ValueError('incomplete twist info')
        elif len(args) == 5:
            self.airfoil = args[0]
            self.chord_len = args[1]
            self.twist_ang = args[2]
            self.twist_center = args[3]
            assert len(self.twist_center) == 3
            self.twist_ref = args[4]
            assert 0.0 <= self.twist_ref <= 1.0
        else:
            raise ValueError('invalid input')

    def __repr__(self):
        ret = 'A {:>6.3f}m-long wing-profile'.format(self.chord_len)
        ret += ' based on {:>16},'.format(self.airfoil)
        ret += ' with{:>6.2f} degrees of twist'.format(self.twist_ang)
        ret += ' referring at the {:>5.1f}% of chord.'.format(self.twist_ref)
        return ret

    @classmethod
    def from_geom_param(cls, *args, **kwargs):
        """
        Construct the profile from geometric descriptions.
        :return: Target profile param representation.
        :rtype: WingProfileParam
        """

        foil = args[0]  # Airfoil object
        length = args[1]  # Chord length
        z_offset = args[2]  # Offset in span-wise direction
        swp_back = args[3]  # Angle of sweep-back
        twist = args[4]  # Angle of twist
        dihedral = args[5]  # Angle of dihedral

        '''Referring coordinates'''
        x_ref = kwargs['x_ref'] if 'x_ref' in kwargs else 0.0
        y_ref = kwargs['y_ref'] if 'y_ref' in kwargs else 0.0
        z_ref = kwargs['z_ref'] if 'z_ref' in kwargs else 0.0

        '''Initial endings'''
        front = np.array([x_ref + z_offset * math.tan(math.radians(swp_back)),
                          y_ref + z_offset * math.tan(math.radians(dihedral)),
                          z_ref + z_offset])
        tail = pnt_pan(front, (length, 0, 0))

        '''Center of twist'''
        twist_ref = kwargs['twist_ref'] if 'twist_ref' in kwargs else 1.0
        assert 0.0 <= twist_ref <= 1.0
        twist_center = share(twist_ref, front, tail)

        return cls(foil, length, twist, twist_center, twist_ref)


class WingProfileList(object):
    def __init__(self, *args):
        """
        Construct several wing profiles in one pass.
        :param args: Geom parameters.
        """

        self.pf = []

        if len(args) == 6:
            '''from intrinsic parameters'''
            # check params
            n = len(args[0])
            for i in range(1, 6):
                assert len(args[i]) == n
            # separate params
            airfoil, z, xf, yf, xt, yt = args
            # construct wing profile
            for k in range(n):
                ending = np.array([[xf[k], yf[k], z[k]], [xt[k], yt[k], z[k]]])
                wp = WingProfile(airfoil[k], ending)
                self.pf.append(wp)
        elif len(args) == 2:
            '''from final parameters'''
            assert len(args[0]) == len(args[1])
            n = len(args[0])
            airfoil, ending = args
            for k in range(n):
                self.pf.append(WingProfile(airfoil[k], ending[k]))
        elif len(args) == 10:
            '''from initial geometric parameters'''
            # check params
            n = len(args[0])
            for i in range(1, 10):
                assert len(args[i]) == n
            # separate params
            airfoil, length, z_off, sweep, twist, dihedral, twist_ref, x_ref, y_ref, z_ref = args
            # construct wing profile
            for k in range(n):
                wpp = WingProfileParam.from_geom_param(airfoil[k], length[k], z_off[k], sweep[k], twist[k], dihedral[k], twist_ref=twist_ref[k], x_ref=x_ref[k], y_ref=y_ref[k], z_ref=z_ref[k])
                wp = WingProfile.from_profile_param(wpp)
                self.pf.append(wp)
        elif len(args) == 5:
            '''from planform'''
            # check params
            n = len(args[1])
            for i in range(2, 5):
                assert len(args[i]) == n
            # separate params
            planform, airfoil, twist_ang, twist_ref, rel_pos = args
            # Construct the wing with given planform
            for i in range(n):
                cu = rel_pos[i]
                cz = planform.z(cu)
                leading = np.array([planform.x_front(cu), planform.y_front(cu), cz])
                trailing = np.array([planform.x_tail(cu), planform.y_tail(cu), cz])
                tst_ref = twist_ref[i]
                tst_center = share(tst_ref, leading, trailing)
                tst_ang = twist_ang[i]
                cur_twist = math.radians(tst_ang)
                actual_len = planform.chord_len(cu) / math.cos(cur_twist)
                wpp = WingProfileParam(airfoil[i], actual_len, tst_ang, tst_center, tst_ref)
                wp = WingProfile.from_profile_param(wpp)
                self.pf.append(wp)
        else:
            raise ValueError('unknown input')

    @property
    def size(self):
        return len(self.pf)

    def add(self, p):
        self.pf.append(p)

    def clear(self):
        self.pf.clear()

    def at(self, idx):
        """
        Refer to an element.
        :param idx: Index of the element.
        :type idx: int
        :return: Profile at given index.
        :rtype: WingProfile
        """

        return self.pf[idx]

    def crv_list_in_nurbs(self):
        ret = []
        for k, elem in enumerate(self.pf):
            cur_crv = elem.crv
            ret.append(cur_crv)
        return ret


class Wing(object):
    def __init__(self, wpl):
        """
        Wing constructed from profiles in span-wise direction.
        :param wpl: List of wing profiles.
        :type wpl: WingProfileList
        """

        self.profile = wpl
        self.surf = self._construct_surf()

    def __repr__(self):
        return "Wing with {} sections".format(self.size)

    @property
    def size(self):
        return self.profile.size

    def at(self, idx):
        """
        Referring wing profile.
        :param idx: Target index.
        :type idx: int
        :return: Target wing profile.
        :rtype: WingProfile
        """

        return self.profile.at(idx)

    @property
    def root(self):
        """
        Profile at root.
        :return: The root profile in WingProfile representation.
        :rtype: WingProfile
        """

        return self.at(0)

    @property
    def tip(self):
        """
        Profile at tip.
        :return: The tip profile in WingProfile representation.
        :rtype: WingProfile
        """

        return self.at(-1)

    def _construct_surf(self):
        return Skinned([self.at(i).crv for i in range(self.size)], 3, 3)

    def _update_surf(self):
        self.surf = self._construct_surf()

    @property
    def leading(self):
        """
        Leading edge of the wing.
        :return: Leading edge in NURBS representation.
        :rtype: Crv
        """

        pts = [self.at(i).front for i in range(self.size)]
        return Spline(pts, method='chord')

    @property
    def tailing_up(self):
        return self.surf.extract('U', 0)

    @property
    def tailing_down(self):
        return self.surf.extract('U', 1)

    def gen_grid(self, *args, **kwargs):
        """
        Generate the multi-block grid for a simplified wing.
        :param args: Containing the geometric description and node distribution of the grid.
        :param kwargs: Extra options on smoothing and spacing.
        :return: The wire-frame, plot3d-grid and fluent-grid(with predefined BC) of the flow field.
        """

        assert len(args) == 5

        a, b, c, d, n = args
        assert len(n) == 8

        iges_model = Model()
        p3d_grid = Plot3D()

        '''
        wsf: Wing surface.
        rt: Wing root curve.
        rl: Wing root chord length.
        tp: Wing tip curve.
        z: Z-Dim position of root, tip and far.
        spn: Wing span(From root to tip).
        far: Wing tip extension on far-field.
        wing_brk: Splitting position on wing surf in U direction to construct separated blocks.
        inlet_brk: Splitting position on inlet surf in U direction to construct separated blocks.
        '''
        wsf = self.surf
        rt = wsf.extract('V', 0)
        rl = self.profile[0].chord_len
        tp = wsf.extract('V', 1)
        z = np.array([self.profile[0].z, self.profile[-1].z, self.profile[-1].z + d])
        spn = z[1] - z[0]
        far = WingProfile.from_geom_param(self.profile[-1].name, z[2], rl, 0, 0, 0, thickness_factor=2).crv
        fsf = RuledSurf(tp, far)
        wing_brk = kwargs['wing_brk_param'] if 'wing_brk_param' in kwargs else (0.45, 0.55)
        inlet_brk = kwargs['inlet_brk'] if 'inlet_brk' in kwargs else (0.3, 0.72)

        assert len(wing_brk) == 2
        assert len(inlet_brk) == 2

        t1 = (rl, b, z[0])
        t2 = (rl, -b, z[0])
        inlet1 = ConicArc(t1, (-1, 0, 0), t2, (1, 0, 0), (rl - a, 0, z[0]))
        t1 = pnt_pan(inlet1.start, (0, 0, spn))
        t2 = pnt_pan(inlet1.end, (0, 0, spn))
        inlet2 = ConicArc(t1, (-1, 0, 0), t2, (1, 0, 0), (rl - a, 0, z[1]))
        t1 = pnt_pan(inlet2.start, (0, 0, d))
        t2 = pnt_pan(inlet2.end, (0, 0, d))
        inlet3 = ConicArc(t1, (-1, 0, 0), t2, (1, 0, 0), (rl - a, 0, z[2]))

        '''Points, lines, curves, surfs'''
        p = np.zeros((36, 3))
        p[0] = inlet1.start
        p[1] = (p[0][0] + c, p[0][1], z[0])
        p[2] = wsf(0, 0)
        p[3] = (p[1][0], p[2][1], z[0])
        p[4] = wsf(1, 0)
        p[5] = (p[3][0], p[4][1], z[0])
        p[6] = inlet1.end
        p[7] = (p[5][0], p[6][1], z[0])
        p[8] = inlet2.start
        p[9] = pnt_pan(p[1], (0, 0, spn))
        p[10] = wsf(0, 1)
        p[11] = (p[9][0], p[10][1], z[1])
        p[12] = wsf(1, 1)
        p[13] = (p[11][0], p[12][1], z[1])
        p[14] = inlet2.end
        p[15] = pnt_pan(p[7], (0, 0, spn))
        p[16] = inlet3.start
        p[17] = pnt_pan(p[9], (0, 0, d))
        p[18] = far.start
        p[19] = (p[17][0], p[18][1], z[2])
        p[20] = far.end
        p[21] = (p[19][0], p[20][1], z[2])
        p[22] = inlet3.end
        p[23] = pnt_pan(p[15], (0, 0, d))
        p[24] = rt(wing_brk[0])
        p[25] = rt(wing_brk[1])
        p[26] = tp(wing_brk[0])
        p[27] = tp(wing_brk[1])
        p[28] = far(wing_brk[0])
        p[29] = far(wing_brk[1])
        p[30] = inlet1(inlet_brk[0])
        p[31] = inlet1(inlet_brk[1])
        p[32] = inlet2(inlet_brk[0])
        p[33] = inlet2(inlet_brk[1])
        p[34] = inlet3(inlet_brk[0])
        p[35] = inlet3(inlet_brk[1])

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
             Line(p[33], p[35])]  # 55

        c = [wsf.extract('U', 0), wsf.extract('U', 1)]  # c0, c1
        c += Crv.split(inlet1, inlet_brk)  # c2, c3, c4
        c += Crv.split(inlet2, inlet_brk)  # c5, c6, c7
        c += Crv.split(inlet3, inlet_brk)  # c8, c9, c10
        c += Crv.split(rt, wing_brk)  # c11, c12, c13
        c += Crv.split(tp, wing_brk)  # c14, c15, c16
        c += Crv.split(far, wing_brk)  # c17, c18, c19
        c += [wsf.extract('U', bp) for bp in wing_brk]  # c20, c21
        c += [fsf.extract('U', bp) for bp in wing_brk]  # c22, c23
        for k in (3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19):
            c[k].reverse()

        ts1 = Surf.split(wsf, wing_brk, [])
        ts2 = Surf.split(fsf, wing_brk, [])
        s = [ts1[0][0], ts1[1][0], ts1[2][0], ts2[0][0], ts2[1][0], ts2[2][0]]

        '''IGES Model'''
        for pnt in p:
            iges_model.add(Entity116(pnt[0], pnt[1], pnt[2]))
        for line in l:
            iges_model.add(Entity110(line.start, line.end))
        for crv in c:
            iges_model.add(crv.to_iges())
        # for surf in s:
        #     iges_model.add(surf.to_iges())
        for elem in self.profile:
            iges_model.add(elem.crv.to_iges())

        '''Construct blocks'''
        blk = [LinearTFI3D.from_edges(l[1], l[26], l[0], l[35], l[5], l[29], l[4], l[37], c[0], l[13], l[12], l[20]),  # BLK0
               LinearTFI3D.from_edges(l[2], l[27], l[1], l[39], l[6], l[30], l[5], l[40], c[1], l[14], l[13], c[0]),  # BLK1
               LinearTFI3D.from_edges(l[36], l[3], l[28], l[2], l[38], l[7], l[31], l[6], c[1], l[21], l[15], l[14]),  # BLK2
               LinearTFI3D.from_edges(l[5], l[29], l[4], l[37], l[9], l[32], l[8], l[41], l[23], l[17], l[16], l[22]),  # BLK3
               LinearTFI3D.from_edges(l[6], l[30], l[5], l[40], l[10], l[33], l[9], l[42], l[24], l[18], l[17], l[23]),  # BLK4
               LinearTFI3D.from_edges(l[38], l[7], l[31], l[6], l[43], l[11], l[34], l[10], l[24], l[25], l[19], l[18])]  # BLK5

        b6_s1 = deepcopy(s[0])
        b6_s2 = LinearTFI2D(c[2], l[20], c[5], l[52])
        b6_s3 = LinearTFI2D(c[0], l[35], l[20], l[37])
        b6_s4 = LinearTFI2D(c[20], l[44], l[52], l[46])
        b6_s5 = LinearTFI2D(l[35], c[11], l[44], c[2])
        b6_s6 = LinearTFI2D(l[37], c[14], l[46], c[5])
        b6_tfi_grid = LinearTFI3D(lambda v, w: b6_s1(v, w), lambda v, w: b6_s2(v, w),
                                  lambda w, u: b6_s3(w, u), lambda w, u: b6_s4(w, u),
                                  lambda u, v: b6_s5(u, v), lambda u, v: b6_s6(u, v))
        blk.append(b6_tfi_grid)

        b7_s1 = LinearTFI2D(l[45], c[21], l[47], l[53])
        b7_s2 = LinearTFI2D(l[44], c[20], l[46], l[52])
        b7_s3 = deepcopy(s[1])
        b7_s3.reverse('U')
        b7_s3.swap()
        b7_s4 = LinearTFI2D(l[53], c[3], l[52], c[6])
        b7_s5 = LinearTFI2D(c[12], l[45], c[3], l[44])
        b7_s6 = LinearTFI2D(c[15], l[47], c[6], l[46])
        b7_tfi_grid = LinearTFI3D(lambda v, w: b7_s1(v, w), lambda v, w: b7_s2(v, w),
                                  lambda w, u: b7_s3(w, u), lambda w, u: b7_s4(w, u),
                                  lambda u, v: b7_s5(u, v), lambda u, v: b7_s6(u, v))
        blk.append(b7_tfi_grid)

        b8_s1 = LinearTFI2D(l[36], c[1], l[38], l[21])
        b8_s2 = LinearTFI2D(l[45], c[21], l[47], l[53])
        b8_s3 = deepcopy(s[2])
        b8_s3.reverse('U')
        b8_s3.swap()
        b8_s4 = LinearTFI2D(l[21], c[4], l[53], c[7])
        b8_s5 = LinearTFI2D(c[13], l[36], c[4], l[45])
        b8_s6 = LinearTFI2D(c[16], l[38], c[7], l[47])
        b8_tfi_grid = LinearTFI3D(lambda v, w: b8_s1(v, w), lambda v, w: b8_s2(v, w),
                                  lambda w, u: b8_s3(w, u), lambda w, u: b8_s4(w, u),
                                  lambda u, v: b8_s5(u, v), lambda u, v: b8_s6(u, v))
        blk.append(b8_tfi_grid)

        b9_s1 = deepcopy(s[3])
        b9_s2 = LinearTFI2D(c[5], l[22], c[8], l[54])
        b9_s3 = LinearTFI2D(l[23], l[37], l[22], l[41])
        b9_s4 = LinearTFI2D(c[22], l[46], l[54], l[48])
        b9_s5 = LinearTFI2D(l[37], c[14], l[46], c[5])
        b9_s6 = LinearTFI2D(l[41], c[17], l[48], c[8])
        b9_tfi_grid = LinearTFI3D(lambda v, w: b9_s1(v, w), lambda v, w: b9_s2(v, w),
                                  lambda w, u: b9_s3(w, u), lambda w, u: b9_s4(w, u),
                                  lambda u, v: b9_s5(u, v), lambda u, v: b9_s6(u, v))
        blk.append(b9_tfi_grid)

        b10_s1 = LinearTFI2D(l[47], c[23], l[49], l[55])
        b10_s2 = LinearTFI2D(l[46], c[22], l[48], l[54])
        b10_s3 = deepcopy(s[4])
        b10_s3.reverse('U')
        b10_s3.swap()
        b10_s4 = LinearTFI2D(l[55], c[6], l[54], c[9])
        b10_s5 = LinearTFI2D(c[15], l[47], c[6], l[46])
        b10_s6 = LinearTFI2D(c[18], l[49], c[9], l[48])
        b10_tfi_grid = LinearTFI3D(lambda v, w: b10_s1(v, w), lambda v, w: b10_s2(v, w),
                                   lambda w, u: b10_s3(w, u), lambda w, u: b10_s4(w, u),
                                   lambda u, v: b10_s5(u, v), lambda u, v: b10_s6(u, v))
        blk.append(b10_tfi_grid)

        b11_s1 = LinearTFI2D(l[38], l[24], l[43], l[25])
        b11_s2 = LinearTFI2D(l[47], c[23], l[49], l[55])
        b11_s3 = deepcopy(s[5])
        b11_s3.reverse('U')
        b11_s3.swap()
        b11_s4 = LinearTFI2D(l[25], c[7], l[55], c[10])
        b11_s5 = LinearTFI2D(c[16], l[38], c[7], l[47])
        b11_s6 = LinearTFI2D(c[19], l[43], c[10], l[49])
        b11_tfi_grid = LinearTFI3D(lambda v, w: b11_s1(v, w), lambda v, w: b11_s2(v, w),
                                   lambda w, u: b11_s3(w, u), lambda w, u: b11_s4(w, u),
                                   lambda u, v: b11_s5(u, v), lambda u, v: b11_s6(u, v))
        blk.append(b11_tfi_grid)

        b12_s1 = deepcopy(s[5])
        b12_s1.reverse('U')
        b12_s2 = deepcopy(s[3])
        b12_s3 = LinearTFI2D(l[24], l[40], l[23], l[42])
        b12_s4 = deepcopy(s[4])
        b12_s4.reverse('U')
        b12_s4.swap()
        b12_s5 = LinearTFI2D(l[40], c[16], c[15], c[14])
        b12_s6 = LinearTFI2D(l[42], c[19], c[18], c[17])
        b12_tfi_grid = LinearTFI3D(lambda v, w: b12_s1(v, w), lambda v, w: b12_s2(v, w),
                                   lambda w, u: b12_s3(w, u), lambda w, u: b12_s4(w, u),
                                   lambda u, v: b12_s5(u, v), lambda u, v: b12_s6(u, v))
        blk.append(b12_tfi_grid)

        '''Node distribution'''
        node = [hyperbolic_tangent(n[0], 8),  # u0
                double_exponential(n[1], 0.5, 1.5, 0.5),  # u1
                uniform(n[2]),  # u2
                single_exponential(n[3], 5),  # u3
                hyperbolic_tangent(n[4], 5),  # u4
                double_exponential(n[5], 0.5, 1.2, 0.5),  # u5
                double_exponential(n[6], 0.5, 1.5, 0.5),  # u6
                uniform(n[7])]  # u7

        blk_node_param = [(3, 0, 2), (3, 7, 2), (0, 3, 2), (3, 0, 4), (3, 7, 4), (0, 3, 4), (0, 1, 2), (5, 0, 2), (6, 0, 2), (0, 1, 4), (5, 0, 4), (6, 0, 4), (7, 6, 4)]
        assert len(blk_node_param) == len(blk)

        '''Calculate grid'''
        # print('Calculating grid...')
        # for i in range(len(blk)):
        #     print('Calculate blk{}...'.format(i))
        #     nu, nv, nw = blk_node_param[i]
        #     blk[i].calc_grid(node[nu], node[nv], node[nw])
        # tfi_grid = [blk[i].grid for i in range(len(blk))]

        '''Smoothing'''
        # print('Smoothing...')
        # for i in (6, 7, 8, 12):
        #     print('Smoothing blk{}...'.format(i))
        #     l3d = Laplace3D(tfi_grid[i])
        #     l3d.smooth()
        #     tfi_grid[i] = np.copy(l3d.grid)

        '''Build Plot3D Output'''
        # for i in range(len(tfi_grid)):
        #     p3d_grid.add(Plot3DBlock.construct_from_array(tfi_grid[i]))

        # return iges_model, p3d_grid
        return iges_model
