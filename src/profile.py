#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
from misc import share
from settings import z_axis_positive
from load_dist import LiftDist
from planform import WingPlanform
from airfoil import Airfoil


def calc_profile_cl(rel_pos, lift_dist, wing_planform, q_inf):
    """
    Approximate calculation of lift coefficient in a finite area.
    :param rel_pos: Relative position list.
    :param lift_dist: Lift distribution. It is assumed to be sorted in ascending order.
    :type lift_dist: LiftDist
    :param wing_planform: Planform of the wing, which provides the b(z).
    :type wing_planform: WingPlanform
    :param q_inf: Dynamic pressure of the free-stream.
    :type q_inf: float
    :return: Approximate Cl.
    """

    n = len(rel_pos)
    lift = np.zeros(n)
    area = np.zeros(n)

    pos = [rel_pos[0]]
    for i in range(1, n):
        mid = 0.5 * (rel_pos[i - 1] + rel_pos[i])
        pos.append(mid)
    pos.append(rel_pos[-1])

    for i in range(n):
        lift[i] = lift_dist.lift_between(pos[i], pos[i + 1])
        area[i] = wing_planform.area_between(pos[i], pos[i + 1])

    return np.array([lift[i] / (q_inf * area[i]) for i in range(n)])


class ProfileSpatialParam(object):
    def __init__(self, chord, twist_ref, twist_center, twist_ang):
        pass

    @classmethod
    def from_leading(cls, leading, chord, incidence):
        pass

    @classmethod
    def from_trailing(cls, trailing, chord, incidence):
        pass

    @classmethod
    def from_2pnt(cls, leading, trailing):
        pass

    @classmethod
    def from_spatial(cls, chord, z_off, twist_ang, dihedral_ang, twist_ref=1.0, x_ref=0, y_ref=0, z_ref=0):
        pass


class Profile(object):
    def __init__(self, airfoil, param):
        """
        Profile of a wing at certain position in span-wise direction.
        :param airfoil: 2D Airfoil.
        :type airfoil: Airfoil
        :param param: Spatial position indicator.
        :type param: ProfileSpatialParam
        """

        self.foil = airfoil
        self.spatial_param = param


class WingProfile(object):
    def __init__(self, foil, ends):
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
