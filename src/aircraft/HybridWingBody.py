import numpy as np
import math
from copy import deepcopy
from matplotlib import pyplot as plt
from src.iges import IGES_Model, IGES_Pnt, IGES_Line
from src.nurbs import ConicArc, LocalCubicInterpolatedCrv, point_inverse
from src.misc import pnt_pan
from src.wing import Wing
from src.aircraft.Baseline import global_origin, z_axis_positive, z_axis_negative
from src.aircraft.Baseline import WingPlanform, VSPlanform, construct_vs_profiles
from src.aircraft.Baseline import HSPlanform, construct_hs_profiles


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

    def __repr__(self):
        return 'Hybrid-Wing-Body Outer Wing Planform'


def construct_hwb_wing_profiles(*args, **kwargs):
    """
    Construct the wing profiles from geometric parameters.
    :param args: Critical description.
    :param kwargs: Optional settings.
    :return: Curves on each profile in a list.
    """

    '''airfoils'''
    foil = args[0]

    '''chord length'''
    chord = np.copy(args[1])

    '''offset in span-wise direction'''
    offset = np.copy(args[2])

    '''angle of sweep-back'''
    sweep_back = np.copy(args[3])  # deg

    '''angle of twist'''
    twist = np.copy(args[4])  # deg

    '''angle of dihedral'''
    dihedral = np.copy(args[5])  # deg

    '''num of profiles'''
    n = len(foil)
    assert n >= 2

    '''other default settings'''
    thickness = np.ones(n, float) if 'thickness' not in kwargs else np.copy(kwargs['thickness'])
    twist_ref = np.ones(n, float) if 'twist_ref' not in kwargs else np.copy(kwargs['twist_ref'])
    dihedral_ref = np.zeros(n, float) if 'dihedral_ref' not in kwargs else np.copy(kwargs['dihedral_ref'])

    '''wing model and profiles'''
    ret = []
    wing = Wing.from_geom_desc(foil, chord, thickness, offset, sweep_back, twist, twist_ref, dihedral, dihedral_ref)
    for elem in wing.profile:
        ret.append(elem.crv)

    '''adjustment'''
    if 'init_origin' in kwargs:
        origin = np.copy(kwargs['init_origin'])
        for i in range(n):
            ret[i].pan(origin)

    if 'incidence' in kwargs:
        incidence_ref = np.copy(kwargs['incidence'][0])
        incidence_ang = kwargs['incidence'][1]
        if not math.isclose(incidence_ang, 0):
            for i in range(n):
                ret[i].rotate(incidence_ref, z_axis_negative, incidence_ang)

    final_delta = [0., 0., 0.]
    if 'downward' in kwargs:
        final_delta[1] -= kwargs['downward']
    if 'forward' in kwargs:
        final_delta[0] -= kwargs['forward']

    for i in range(n):
        ret[i].pan(final_delta)

    return ret


def construct_hwb_frame():
    """
    Generate the Hybrid-Wing-Body model parametrically.
    :return: None.
    """

    '''Fuselage'''
    fuselage_len = 28.0
    span = 42.0
    span2 = span / 2

    fuselage_height = 3.5
    fuselage_width = 5.2

    nose_len = 5.6
    tail_len = 7.8
    body_len = fuselage_len - (tail_len + nose_len)
    print('\nLength of fuselage components:')
    print('Nose: {:.3f}\nBody: {:.3f}\nTail: {:.3f}'.format(nose_len, body_len, tail_len))

    theta_fc = math.radians(15)

    nose_front_downward = 0.8
    nose_front_height = 0.2
    nose_front_width = 0.3
    nose_front_center = pnt_pan(global_origin, (-nose_len, -nose_front_downward, 0))
    nose_front_up = pnt_pan(nose_front_center, (0, nose_front_height / 2, 0))
    nose_front_down = pnt_pan(nose_front_center, (0, -nose_front_height / 2, 0))
    nose_front_mid = pnt_pan(nose_front_center, (0, 0, nose_front_width / 2))
    nose_front_crv = ConicArc(nose_front_up, z_axis_positive, nose_front_down, z_axis_negative, nose_front_mid)

    body_front_up = np.array([0, fuselage_height / 2, 0])
    body_front_down = np.array([0, -fuselage_height / 2, 0])
    body_front_mid = np.array([0, 0, fuselage_width / 2])
    body_front_crv = ConicArc(body_front_up, z_axis_positive, body_front_down, z_axis_negative, body_front_mid)
    body_back_crv = deepcopy(body_front_crv)
    body_front2back_dir = np.array([body_len, 0, 0])
    body_back_crv.pan(body_front2back_dir)
    body_back_up = pnt_pan(body_front_up, body_front2back_dir)
    body_back_down = pnt_pan(body_front_down, body_front2back_dir)
    body_back_mid = pnt_pan(body_front_mid, body_front2back_dir)

    tail_back_height = 0.6
    tail_back_width = 0.6
    tail_back_down = pnt_pan(body_back_down, (tail_len, tail_len * math.tan(theta_fc), 0))
    tail_back_up = pnt_pan(tail_back_down, (0, tail_back_height, 0))
    tail_back_center = 0.5 * (tail_back_up + tail_back_down)
    tail_back_mid = pnt_pan(tail_back_center, (0, 0, tail_back_width / 2))
    tail_back_crv = ConicArc(tail_back_up, z_axis_positive, tail_back_down, z_axis_negative, tail_back_mid)

    '''Wing'''
    fusion_width = 0.4
    wing_spn2 = span2 - fuselage_width / 2 - fusion_width
    wing_root_len = 15.6
    wing_tip_len = 1.5
    wing_leading_inner_delta = (0.8, 1.5)
    wing_leading_middle_delta = (1.25, 1.7)
    wing_leading_outer_sweep = 28
    wing_trailing_inner_delta = (0.7, -2.3)
    wing_trailing_outer_spn = 13.5
    wing_trailing_outer_sweep = 12.5
    wing_planform = HWBWingPlanform(wing_root_len, wing_tip_len, wing_spn2,
                                    wing_leading_inner_delta, wing_leading_middle_delta, wing_leading_outer_sweep,
                                    wing_trailing_inner_delta, wing_trailing_outer_spn, wing_trailing_outer_sweep)

    wing_u = np.array([0.00, 3.33, 6.11, 12.22, 19.44, 28.33, 46.11, 82.22, 100.00]) / 100
    wing_z = wing_u * wing_spn2
    wing_n = len(wing_u)
    wing_ref_origin = (body_len - wing_root_len, 0, span2 - wing_spn2)
    wing_incidence_ref = (body_len + tail_len, 0, fuselage_width / 2)
    wing_incidence_ang = 1.5
    wing_forward_marching = 0.2
    wing_downward_marching = 0.4
    wing_inner_profile_num = 3
    wing_middle_profile_num = 2
    wing_outer_profile_num = wing_n - (wing_inner_profile_num + wing_middle_profile_num)
    wing_foil = ['SC(2)-0710', 'SC(2)-0710', 'SC(2)-0710',
                 'SC(2)-0610', 'SC(2)-0610',
                 'SC(2)-0610', 'SC(2)-0610', 'SC(2)-0410', 'SC(2)-0410']
    wing_tc = np.array([0.01 * float(_f[-2:]) for _f in wing_foil])
    wing_chord = [wing_planform.chord_len(u) for u in wing_u]
    wing_height = [wing_chord[i] * wing_tc[i] for i in range(wing_n)]
    wing_swp = [math.degrees(math.atan2(wing_planform.x_front(u), wing_planform.z(u))) for u in wing_u]
    wing_cl = [0.11, 0.13, 0.15,
               0.70, 0.60,
               0.52, 0.45, 0.30, 0.00]
    wing_twist = [0, -0.1, -0.2,
                  -0.45, -0.95,
                  -2.2, -2.544, -1.765, -3.220]
    # wing_dihedral = [math.atan2((wing_height[i] - wing_height[0]) / 2, wing_z[i]) for i in range(wing_n)]
    wing_dihedral = np.array([2.5, 2.5, 2.5, 3.0, 3.0, 3.5, 3.5, 3.0, 3.0])
    wing_profile = construct_hwb_wing_profiles(wing_foil, wing_chord, wing_z, wing_swp, wing_twist, wing_dihedral,
                                               init_origin=wing_ref_origin,
                                               incidence=[wing_incidence_ref, wing_incidence_ang],
                                               forward=wing_forward_marching, downward=wing_downward_marching)

    print('\nWing:')
    print('Root chord: {:.3f}'.format(wing_root_len))
    print('Half span: {:.3f}'.format(wing_spn2))
    print('Area: {:.2f}'.format(wing_planform.area))
    print('MAC: {:.3f}'.format(wing_planform.mean_aerodynamic_chord))
    print('Dihedral: {}'.format(wing_dihedral))

    '''Vertical-Stabilizer'''
    vs_root_chord = 5.5
    vs_tip_chord = 3.5
    vs_spn2 = 4.2
    vs_leading_swp = 45
    vs_planform = VSPlanform(vs_root_chord, vs_tip_chord, vs_spn2, vs_leading_swp)
    vs_u = np.array([0.00, 6.67, 13.33, 27.50, 55.00, 85.00, 100.00]) / 100
    vs_n = len(vs_u)
    vs_z = np.array([vs_planform.z(_u) for _u in vs_u])

    print('\nVertical-Stabilizer:')
    print('Area: {:.2f}'.format(vs_planform.area / 2))
    print('MAC: {:.3f}'.format(vs_planform.mean_aerodynamic_chord))

    vs_foil = ['NACA0010'] * vs_n
    vs_cl = [vs_planform.chord_len(vs_u[i]) for i in range(vs_n)]
    vs_swp = [math.degrees(math.atan2(vs_planform.x_front(vs_u[i]), vs_z[i])) for i in range(vs_n)]

    vs_delta_tail_back = 0.8
    vs_delta_body_up = 0
    vs_delta_x = body_len + tail_len - (vs_root_chord + vs_delta_tail_back)
    vs_delta_y = fuselage_height / 2 + vs_delta_body_up
    vs_origin = (vs_delta_x, vs_delta_y, 0)
    vs_profile = construct_vs_profiles(vs_foil, vs_cl, vs_z, vs_swp, origin=vs_origin)

    '''Horizontal-Stabilizer'''
    hs_root_chord = 3.5
    hs_tip_chord = 1.7
    hs_spn2 = 6.5
    hs_leading_swp = 25
    hs_planform = HSPlanform(hs_root_chord, hs_tip_chord, hs_spn2, hs_leading_swp)
    hs_u = np.array([0.00, 50.00, 100.00]) / 100
    hs_z = hs_u * hs_spn2
    hs_n = len(hs_u)

    print('\nHorizontal-Stabilizer:')
    print('Area: {:.2f}'.format(hs_planform.area))
    print('MAC: {:.3f}'.format(hs_planform.mean_aerodynamic_chord))

    hs_foil = ['NACA0008'] * hs_n
    hs_cl = [hs_planform.chord_len(hs_u[i]) for i in range(hs_n)]
    hs_swp = [math.degrees(math.atan2(hs_planform.x_front(hs_u[i]), hs_z[i])) for i in range(hs_n)]

    hs_forward_marching_ratio = 0
    hs_delta_x = vs_origin[0] + vs_planform.x_front(1) - hs_forward_marching_ratio * vs_tip_chord
    hs_delta_y = vs_origin[1] + vs_spn2
    hs_origin = (hs_delta_x, hs_delta_y, 0)
    hs_profile = construct_hs_profiles(hs_foil, hs_cl, hs_z, hs_swp, origin=hs_origin)

    '''Graphic representation'''
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    wing_planform.pic(ax1, u=wing_u)
    ax2 = fig.add_subplot(222)
    vs_planform.pic(ax2, u=vs_u, direction='vertical')
    # ax3 = fig.add_subplot(223, sharex=ax1)
    # ax3.plot(wing_z, wing_cl, label='Cl')
    # ax3.legend()
    ax4 = fig.add_subplot(224)
    hs_planform.pic(ax4, u=hs_u)
    fig.tight_layout()
    plt.show()

    '''CAD model'''
    model = IGES_Model()
    model.add(IGES_Pnt(body_front_up))
    model.add(IGES_Pnt(body_front_down))
    model.add(IGES_Pnt(body_front_mid))
    model.add(IGES_Line(body_front_up, body_back_up))
    model.add(IGES_Line(body_front_down, body_back_down))
    model.add(IGES_Line(body_front_mid, body_back_mid))
    model.add(body_front_crv.to_iges())
    model.add(body_back_crv.to_iges())
    model.add(tail_back_crv.to_iges())
    model.add(IGES_Line(tail_back_center, tail_back_up))
    model.add(IGES_Line(tail_back_center, tail_back_down))
    model.add(IGES_Line(tail_back_center, tail_back_mid))
    model.add(nose_front_crv.to_iges())
    model.add(IGES_Line(nose_front_center, nose_front_up))
    model.add(IGES_Line(nose_front_center, nose_front_down))
    model.add(IGES_Line(nose_front_center, nose_front_mid))
    for _c in wing_profile:
        model.add(_c.to_iges())
    for _c in vs_profile:
        model.add(_c.to_iges())
    for _c in hs_profile:
        model.add(_c.to_iges())
    model.save('HWB.igs')


if __name__ == '__main__':
    construct_hwb_frame()
