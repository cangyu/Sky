from matplotlib import pyplot as plt
from abc import abstractmethod, ABCMeta
from copy import deepcopy
import math
import numpy as np
from numpy.linalg import norm
from scipy.integrate import romberg
from grid import LinearTFI2D, LinearTFI3D, Laplace2D, ThomasMiddlecoff2D
from grid import Plot3D, Plot3DBlock
from grid import hyperbolic_tangent, uniform, single_exponential, double_exponential
from iges import Model, Entity116, Entity110
from misc import pnt_dist, read_airfoil_pts, pnt_pan, share
from nurbs import Crv, Line, Spline, ConicArc
from nurbs import GlobalInterpolatedCrv, LocalCubicInterpolatedCrv
from nurbs import Surf, Skinned, RuledSurf, point_inverse
from settings import AIRFOIL_LIST

global_origin = np.array([0., 0., 0.])
x_axis_positive = np.array([1., 0, 0])
x_axis_negative = np.array([-1., 0, 0])
y_axis_positive = np.array([0, 1., 0])
y_axis_negative = np.array([0, -1., 0])
z_axis_positive = np.array([0, 0, 1.])
z_axis_negative = np.array([0, 0, -1.])


class EllipticLiftDist(object):
    def __init__(self, w, spn, rho_inf, vel_inf):
        """
        Elliptic Lift distribution in span-wise direction.
        This is the ideal distribution for reducing induced drag.
        :param w: Weight of the aircraft in cruise, in tons.
        :type w: float
        :param spn: Span of the aircraft, in meters.
        :type spn: float
        :param rho_inf: Density of the free-stream.
        :type rho_inf: float
        :param vel_inf: Velocity of the free-stream.
        :type vel_inf: float
        """

        self.payload = w * 1000 * 9.8 / 2
        self.span2 = spn / 2
        self.rho = rho_inf
        self.v = vel_inf
        self.p_inf = 0.5 * self.rho * self.v ** 2
        self.root_lift = self.payload / (0.25 * math.pi * self.span2)

    def lift_at(self, rel_pos):
        return self.root_lift * math.sqrt(1 - rel_pos ** 2)

    def velocity_circulation_at(self, rel_pos):
        return self.lift_at(rel_pos) / (self.rho * self.v)

    def cl_at(self, rel_pos, chord_len):
        return self.lift_at(rel_pos) / (self.p_inf * chord_len)


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

        return 2 * math.fabs(romberg(self.chord_len, 0, 1)) * self.z(1)

    @property
    def mean_aerodynamic_chord(self):
        """
        Get the mean aerodynamic chord length of the wing.
        :return: The MAC.
        :rtype: float
        """

        return math.fabs(romberg(lambda u: self.chord_len(u) ** 2, 0, 1)) / (0.5 * self.area) * self.z(1)

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
    def __repr__(self):
        pass

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


class HWBWingPlanform(WingPlanform):
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

    def __repr__(self):
        return 'Hybrid-Wing-Body Outer Wing Planform'


class Airfoil(object):
    def __init__(self, foil):
        """
        2D Airfoil, with chord length equals to 1.
        :param foil: Airfoil name(Capital Case).
        :type foil: str
        """

        if foil not in AIRFOIL_LIST:
            raise FileNotFoundError('Airfoil \'{}\' not included at present.'.format(foil))

        self.name = foil
        self.pts = read_airfoil_pts(foil)

    def __repr__(self):
        return '{} with {} points'.format(self.name, self.pnt_num)

    @property
    def z(self):
        return self.pts[0][2]

    @property
    def pnt_num(self):
        return len(self.pts)

    @property
    def crv(self):
        # return Spline(self.pts)
        return GlobalInterpolatedCrv(self.pts, 3)

    @property
    def tail_up(self):
        return self.pts[0]

    @property
    def tail_down(self):
        return self.pts[-1]

    @property
    def tail(self):
        return (self.tail_up + self.tail_down) / 2

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
    def chord_len(self):
        return pnt_dist(self.front, self.tail)

    @property
    def is_blunt(self):
        return not math.isclose(norm(self.tail_up - self.tail_down), 0)

    def to_blunt(self):
        pfx = self.front[0]
        r0 = self.chord_len
        self.pts = self.pts[1:-1]
        r1 = self.chord_len
        ratio = r0 / r1
        for k in range(len(self.pts)):
            self.pts[k][0] = pfx + (self.pts[k][0] - pfx) * ratio

    def save(self, fn):
        """
        Save all coordinates into file.
        :param fn: File name.
        :type fn: str
        :return: None.
        """

        f_out = open(fn, 'w')
        for p in self.pts:
            f_out.write('{:10.6f}\t{:10.6f}\t{:10.6f}\n'.format(p[0], p[1], p[2]))
        f_out.close()

    def plot(self, ax):
        (px, py, pz) = zip(*self.pts)
        ax.plot(px, py)
        ax.set_aspect('equal')

    def curvature_at(self, rel_pos):
        """
        Calculate the curvature at given position.
        :param rel_pos: Relative position.
        :type rel_pos: float
        :return: Curvature.
        :rtype: float
        """

        return self.curve.curvature(rel_pos)

    def gen_grid(self, *args, **kwargs):
        """
        Generate grid for 2D airfoil or wing profile.
        :param args: Containing the geometric description and node distribution of the grid.
        :param kwargs: Extra options on smoothing and spacing.
        :return: The wire-frame, plot3d-grid and fluent-grid(with predefined BC) of the flow field.
        """

        '''
        a: Width of the front part of the flow field.
        b: Semi-height of the flow field.
        c: Width of the rear part of the flow field.
        n0: Num of points along airfoil.
        n1: Num of points along vertical direction.
        n2: Num of points along horizontal direction in rear field.
        n3: Num of Points on the trailing edge.
        '''
        assert len(args) == 7
        a, b, c, n0, n1, n2, n3 = args

        wire_frame = Model()
        p3d_grid = Plot3D()

        '''Flow-field geometries'''
        pts = np.empty((8, 3), float)
        pts[0] = self.tail_up
        pts[1] = self.tail_down
        pts[2] = np.array([0, b, self.z])
        pts[3] = np.array([0, -b, self.z])
        pts[4] = np.array([c, pts[0][1], self.z])
        pts[5] = np.array([c, pts[1][1], self.z])
        pts[6] = np.array([c, b, self.z])
        pts[7] = np.array([c, -b, self.z])

        crv = [self.crv,  # c0
               ConicArc(pts[2], (-1, 0, 0), pts[3], (1, 0, 0), (-a, 0, 0)),  # c1
               Line(pts[0], pts[2]),  # c2
               Line(pts[1], pts[3]),  # c3
               Line(pts[4], pts[6]),  # c4
               Line(pts[5], pts[7]),  # c5
               Line(pts[0], pts[4]),  # c6
               Line(pts[1], pts[5]),  # c7
               Line(pts[2], pts[6]),  # c8
               Line(pts[3], pts[7]),  # c9
               Line(pts[0], pts[1]),  # c10
               Line(pts[4], pts[5])]  # c11

        '''Construct wire-frame'''
        for p in pts:
            wire_frame.add(Entity116(p[0], p[1], p[2]))
        for c in crv:
            wire_frame.add(c.to_iges())

        '''Knot distribution'''
        u = [double_exponential(n0, 0.5, -1.5, 0.5),  # c0, c1
             hyperbolic_tangent(n1, 2),  # c2, c3, c4, c5
             single_exponential(n2, 3),  # c6, c7, c8, c9
             uniform(n3)]  # c10, c11

        '''Structured grid blocks'''
        leading_blk = LinearTFI2D(crv[2], crv[0], crv[3], crv[1])
        tailing_up_blk = LinearTFI2D(crv[6], crv[2], crv[8], crv[4])
        tailing_down_blk = LinearTFI2D(crv[3], crv[7], crv[5], crv[9])
        rear_blk = LinearTFI2D(crv[10], crv[6], crv[11], crv[7])

        '''Construct Plot3D grid for basic checking'''
        leading_blk.calc_grid(u[1], u[0])
        leading_tfi_grid = leading_blk.grid
        leading_smooth_ok = False
        if 'leading_smooth' in kwargs:
            smooth = kwargs['leading_smooth']
            if smooth in ('Laplace', 'laplace'):
                leading_grid_laplace = Laplace2D(leading_tfi_grid)
                leading_grid_laplace.smooth()
                p3d_grid.add(Plot3DBlock.construct_from_array(leading_grid_laplace.grid))
                leading_smooth_ok = True
            if smooth in ('TM', 'tm', 'Thomas-Middlecoff', 'thomas-middlecoff'):
                leading_grid_tm = ThomasMiddlecoff2D(leading_tfi_grid)
                leading_grid_tm.smooth()
                p3d_grid.add(Plot3DBlock.construct_from_array(leading_grid_tm.grid))
                leading_smooth_ok = True
        if not leading_smooth_ok:
            p3d_grid.add(Plot3DBlock.construct_from_array(leading_tfi_grid))

        tailing_up_blk.calc_grid(u[2], u[1])
        tailing_up_grid = tailing_up_blk.grid
        p3d_grid.add(Plot3DBlock.construct_from_array(tailing_up_grid))

        tailing_down_blk.calc_grid(u[1], u[2])
        tailing_down_grid = tailing_down_blk.grid
        p3d_grid.add(Plot3DBlock.construct_from_array(tailing_down_grid))

        if self.is_blunt:
            rear_blk.calc_grid(u[3], u[2])
            rear_grid = rear_blk.grid
            p3d_grid.add(Plot3DBlock.construct_from_array(rear_grid))

        return wire_frame, p3d_grid


class WingProfileParam(object):
    def __init__(self, *args, **kwargs):
        """
        Intrinsic descriptions for profile in span-wise direction.
        :param args: Elementary parameters.
        :param kwargs: Optional parameters.
        """

        # Name of the airfoil, input in str
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

        # Relative thickness, input in percentage [0, 100]
        self.tc = kwargs['thickness'] if 'thickness' in kwargs else float(self.airfoil[-2:])
        self.tc *= 0.01
        assert 0.0 <= self.tc <= 1.0

        # Maximum height
        self.height = self.tc * self.chord_len

        # Aerodynamic properties
        self.cl = 0.0 if 'cl' not in kwargs else kwargs['cl']
        self.cd = 0.0 if 'cd' not in kwargs else kwargs['cd']
        self.cm = 0.0 if 'cm' not in kwargs else kwargs['cm']

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

        foil = args[0]  # Name of the airfoil
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


class WingProfile(Airfoil):
    def __init__(self, foil, ends):
        """
        3D profile at certain position.
        :param foil: Airfoil name.
        :type foil: str
        :param ends: Starting and ending points of the profile.
        """

        super(WingProfile, self).__init__(foil)

        '''Inspect endings'''
        self.ending = np.copy(ends)
        if not math.isclose(ends[0][2], ends[1][2]):
            raise AssertionError("Inconsistent ending coordinates in Z direction!")

        cl = self.chord_len
        if math.isclose(cl, 0):
            raise ZeroDivisionError("Invalid ending coordinates in XY direction!")

        rotation = complex((ends[1][0] - ends[0][0]) / cl, (ends[1][1] - ends[0][1]) / cl)

        '''Build profile'''
        if not self.is_blunt:
            self.to_blunt()
        for i in range(self.pnt_num):
            '''Stretch, Z offset and Thickness'''
            self.pts[i][0] *= cl
            self.pts[i][1] *= cl
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
        return [elem.crv for elem in self.pf]


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
