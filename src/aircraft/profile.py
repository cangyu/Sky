#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
from settings import z_axis_negative
from aircraft.load_dist import LiftDist
from aircraft.planform import WingPlanform
from cad.rotation import pnt_rotate
from misc import pnt_pan, pnt_dist, share
from cad.nurbs import Spline


def calc_profile_cl(rel_pos, lift_dist, wing_planform):
    """
    Approximate calculation of lift coefficient in a finite area.
    :param rel_pos: Relative position list. It is assumed to be sorted in ascending order.
    :param lift_dist: Lift distribution.
    :type lift_dist: LiftDist
    :param wing_planform: Planform of the wing, which provides the b(z).
    :type wing_planform: WingPlanform
    :return: Approximate Cl.
    """

    n = len(rel_pos)
    lift = np.zeros(n)
    area = np.zeros(n)
    q_inf = lift_dist.q_inf

    pos = [rel_pos[0]]
    for i in range(1, n):
        mid = 0.5 * (rel_pos[i - 1] + rel_pos[i])
        pos.append(mid)
    pos.append(rel_pos[-1])

    for i in range(n):
        lift[i] = lift_dist.lift_between(pos[i], pos[i + 1])
        area[i] = wing_planform.area_between(pos[i], pos[i + 1])

    return np.array([lift[i] / (q_inf * area[i]) for i in range(n)])


def pic_profile_gamma_cl(vc_ax, cl_ax, lift_dist: LiftDist, planform: WingPlanform, n=1000, factor=1.1):
    x_sp = np.linspace(0, planform.half_span, n)
    u_sp = np.linspace(0, 1, n)

    vc = [lift_dist.gamma_at(u) for u in u_sp]
    cl3 = calc_profile_cl(u_sp, lift_dist, planform)
    swp25 = np.array([planform.swp_025(u) for u in u_sp])
    cl2 = np.array([factor * cl3[i] / math.cos(math.radians(swp25[i])) ** 2 for i in range(n)])

    vc_ax.plot(x_sp, vc, 'r')
    vc_ax.set_xlabel('Span-wise position')
    vc_ax.set_ylabel('Velocity Circulation/({})'.format(r'$m^2 \cdot s^{-1}$'))

    cl_ax.plot(x_sp, cl3, label='Cl in 3D')
    cl_ax.plot(x_sp, cl2, label='Cl in 2D')
    cl_ax.set_ylabel('Lift coefficient')
    cl_ax.legend()


class ProfileSpatialParam(object):
    def __init__(self, chord, twist_ang, twist_center, twist_ref):
        self.chord = chord
        self.twist_ang = twist_ang
        self.twist_center = np.copy(twist_center)
        self.twist_ref = twist_ref
        la = twist_ref * chord
        lb = (1.0 - twist_ref) * chord
        self.leading = pnt_pan(self.twist_center, (-la, 0, 0))
        self.trailing = pnt_pan(self.twist_center, (lb, 0, 0))
        self.leading = pnt_rotate(self.twist_center, z_axis_negative, twist_ang, self.leading)
        self.trailing = pnt_rotate(self.twist_center, z_axis_negative, twist_ang, self.trailing)

    @property
    def z(self):
        return self.twist_center[2]

    @property
    def incidence(self):
        return self.twist_ang

    @classmethod
    def from_leading(cls, leading, chord, incidence):
        return cls(chord, incidence, leading, 0.0)

    @classmethod
    def from_trailing(cls, trailing, chord, incidence):
        return cls(chord, incidence, trailing, 1.0)

    @classmethod
    def from_2pnt(cls, leading, trailing, twist_ref=1.0):
        leading = np.copy(leading)
        trailing = np.copy(trailing)
        assert len(leading) == len(trailing) == 3
        assert leading[2] == trailing[2]
        chord = pnt_dist(leading, trailing)
        delta = leading - trailing
        twist_ang = math.degrees(math.atan2(delta[1], -delta[0]))
        twist_center = share(twist_ref, leading, trailing)
        return cls(chord, twist_ang, twist_center, twist_ref)

    @classmethod
    def from_spatial_desc(cls, planar_chord, z_off, swp, twist, dihedral, twist_ref=1.0, x_ref=0, y_ref=0, z_ref=0):
        """
        There's a little bit difference from constructors above as the given chord is not the final chord,
        and this is due to the way we view these given parameters. Here the parameters are basically viewed
        as planar, and the twist effect is not considered at initial setup.
        :param planar_chord: Initial/Planar chord length.
        :type planar_chord: float
        :param z_off: Offset in Z direction.
        :type z_off: float
        :param swp: Angle of sweep-back.
        :type swp: float
        :param twist: Angle of twist/incidence.
        :type twist: float
        :param dihedral: Angle of dihedral.
        :type dihedral: float
        :param twist_ref: Relative position of the center of twist along the chord.
        :type twist_ref: float
        :param x_ref: Reference ordinate in X direction.
        :type x_ref: float
        :param y_ref:Reference ordinate in X direction.
        :type y_ref: float
        :param z_ref:Reference ordinate in X direction.
        :type z_ref: float
        :return: Description of the spatial configuration.
        :rtype: ProfileSpatialParam
        """

        lx = x_ref + z_off * math.tan(math.radians(swp))
        ly = y_ref + z_off * math.tan(math.radians(dihedral))
        lz = z_ref + z_off
        leading = np.array([lx, ly, lz])
        trailing = pnt_pan(leading, (planar_chord, 0, 0))
        twist_center = share(twist_ref, leading, trailing)
        chord = planar_chord / math.cos(math.radians(twist))
        return cls(chord, twist, twist_center, twist_ref)


def _build_profile_pts(foil, desc):
    """
    Compute the real profile coordinates from airfoil.
    :param foil: Ref airfoil.
    :type foil: Airfoil
    :param desc: Spatial position and pose description.
    :type desc: ProfileSpatialParam
    :return: 3D Points describing the profile in XFOIL format.
    """

    '''Scale airfoil to real profile'''
    chord = desc.chord
    scaled = np.array([(p[0] * chord, p[1] * chord, desc.z) for p in foil.pts])

    '''Move to trailing'''
    cur_trailing = (scaled[0] + scaled[-1]) / 2
    pan_dir = desc.trailing - cur_trailing
    paned = np.array([pnt_pan(p, pan_dir) for p in scaled])

    '''Rotate around trailing'''
    return pnt_rotate(desc.trailing, z_axis_negative, desc.twist_ang, paned)


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
        self.pts = _build_profile_pts(airfoil, param)

    def generate_nurbs_crv(self):
        return Spline(self.pts)


class ProfileList(object):
    def __init__(self):
        """
        Construct several wing profiles in one pass.
        """

        self.pf = []

    @property
    def size(self):
        return len(self.pf)

    def add(self, p):
        """
        Append profile to the container.
        :param p: Target profile.
        :type p: Profile
        :return: None.
        """

        self.pf.append(p)

    def clear(self):
        self.pf.clear()

    def at(self, idx):
        """
        Refer to an element.
        :param idx: Index of the element.
        :type idx: int
        :return: Profile at given index.
        :rtype: Profile
        """

        return self.pf[idx]

    def generate_nurbs_crv_list(self):
        return [elem.generate_nurbs_crv() for elem in self.pf]

    @classmethod
    def from_intrinsic_param(cls, airfoil, z, xf, yf, xt, yt):
        assert len(airfoil) == len(z) == len(xf) == len(yf) == len(xt) == len(yt)
        n = len(airfoil)
        ret = cls()
        for k in range(n):
            leading = (xf[k], yf[k], z[k])
            trailing = (xt[k], yt[k], z[k])
            param = ProfileSpatialParam.from_2pnt(leading, trailing)
            wp = Profile(airfoil[k], param)
            ret.add(wp)
        return ret

    @classmethod
    def from_final_param(cls, airfoil, ending):
        assert len(airfoil) == len(ending)
        n = len(airfoil)
        ret = cls()
        for k in range(n):
            param = ProfileSpatialParam.from_2pnt(ending[0], ending[1])
            wp = Profile(airfoil[k], param)
            ret.add(wp)
        return ret

    @classmethod
    def from_initial_geometric_param(cls, airfoil, length, z_off, sweep, twist, dihedral, twist_ref, x_ref, y_ref, z_ref):
        assert len(airfoil) == len(length) == len(z_off) == len(sweep) == len(twist) == len(dihedral) == len(twist_ref) == len(x_ref) == len(y_ref) == len(z_ref)
        n = len(airfoil)
        ret = cls()
        for k in range(n):
            param = ProfileSpatialParam.from_spatial_desc(length[k], z_off[k], sweep[k], twist[k], dihedral[k], twist_ref[k], x_ref[k], y_ref[k], z_ref[k])
            wp = Profile(airfoil[k], param)
            ret.add(wp)
        return ret

    @classmethod
    def from_planform(cls, planform, airfoil, twist_ang, twist_ref, rel_pos):
        """
        Build wing from planform description.
        :param planform: Planform description.
        :type planform: WingPlanform
        :param airfoil: Airfoil on each profile.
        :param twist_ang: Twist angle of each profile.
        :param twist_ref: Relative twist position on each chord.
        :param rel_pos: Span-wise position for each chord along the planform.
        :return: Target profile list.
        """

        n = len(airfoil)
        assert n == len(twist_ang)
        assert n == len(twist_ref)
        assert n == len(rel_pos)

        ret = cls()
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
            param = ProfileSpatialParam(actual_len, tst_ang, tst_center, tst_ref)
            wp = Profile(airfoil[i], param)
            ret.add(wp)
        return ret

    @classmethod
    def from_planform_with_dihedral(cls, planform, rel_pos, airfoil, twist_ang, twist_ref, dihedral_offset):
        """
        Build wing from the planform, twists and dihedral offsets.
        :param planform: Planform description.
        :type planform: WingPlanform
        :param rel_pos: Span-wise position for each chord along the planform.
        :param airfoil: Airfoil on each profile.
        :param twist_ang: Twist angle of each profile.
        :param twist_ref: Relative twist position on each chord.
        :param dihedral_offset: Offset in Y-Direction relative to the planform surface for each profile.
        :return: Target profile list.
        """

        n = len(airfoil)
        assert n == len(twist_ang)
        assert n == len(twist_ref)
        assert n == len(rel_pos)
        assert n == len(dihedral_offset)

        ret = cls()
        for i in range(n):
            cu = rel_pos[i]
            cz = planform.z(cu)
            leading = pnt_pan([planform.x_front(cu), planform.y_front(cu), cz], [0, dihedral_offset[i], 0])
            trailing = pnt_pan([planform.x_tail(cu), planform.y_tail(cu), cz], [0, dihedral_offset[i], 0])
            tst_ref = twist_ref[i]
            tst_center = share(tst_ref, leading, trailing)
            tst_ang = twist_ang[i]
            cur_twist = math.radians(tst_ang)
            actual_len = planform.chord_len(cu) / math.cos(cur_twist)
            param = ProfileSpatialParam(actual_len, tst_ang, tst_center, tst_ref)
            wp = Profile(airfoil[i], param)
            ret.add(wp)
        return ret
