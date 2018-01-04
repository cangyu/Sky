import numpy as np
import math
import platform
from copy import deepcopy
from matplotlib import pyplot as plt
from iges import IGES_Model, IGES_Pnt, IGES_Line
from nurbs import ConicArc
from misc import pnt_pan, share
from wing import Wing, WingProfileParam, EllipticLiftDist, HWBWingPlanform
from wing import global_origin, z_axis_positive, z_axis_negative
from aircraft import VSPlanform, construct_vs_profiles, HSPlanform, construct_hs_profiles


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

    payload = 120
    ma = 0.75
    a = 299.5
    rho = 0.4135
    velocity = ma * a

    '''Fuselage'''
    fuselage_len = 28.0
    span = 42.0
    span2 = span / 2

    fuselage_height = 3.5
    fuselage_width = 4.0

    nose_len = 4.5
    tail_len = 8.5
    body_len = fuselage_len - (tail_len + nose_len)
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
    tail_back_center = np.copy(share(0.5, tail_back_up, tail_back_down))
    tail_back_mid = pnt_pan(tail_back_center, (0, 0, tail_back_width / 2))
    tail_back_crv = ConicArc(tail_back_up, z_axis_positive, tail_back_down, z_axis_negative, tail_back_mid)

    '''Wing'''
    fusion_width = 0.3
    wing_spn2 = span2 - fuselage_width / 2 - fusion_width
    wing_root_len = 18.4
    wing_tip_len = 1.8
    wing_leading_inner_delta = (1.0, 1.55)
    wing_leading_middle_delta = (1.3, 1.75)
    wing_leading_outer_sweep = 28
    wing_trailing_inner_delta = (0.7, -2.3)
    wing_trailing_outer_spn = 13.5
    wing_trailing_outer_sweep = 12
    wing_planform = HWBWingPlanform(wing_root_len, wing_tip_len, wing_spn2,
                                    wing_leading_inner_delta, wing_leading_middle_delta, wing_leading_outer_sweep,
                                    wing_trailing_inner_delta, wing_trailing_outer_spn, wing_trailing_outer_sweep)

    wing_u = np.array([0.0000, 0.0400, 0.0985, 0.1793, 0.2600,
                       0.3850, 0.5100, 0.6200, 0.7400, 0.8500, 0.9300, 0.9800, 1.0000])
    wing_n = len(wing_u)
    wing_chord = np.array([wing_planform.chord_len(_u) for _u in wing_u])
    wing_ref_origin = (body_len - wing_root_len, 0, span2 - wing_spn2)

    wing_swp = [math.degrees(math.atan2(wing_planform.x_front(u), wing_planform.z(u))) for u in wing_u]
    wing_dihedral = np.array([2.5, 2.5, 3.0, 3.0, 3.5, 3.5, 3.5, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0])
    wing_z = wing_ref_origin[2] + wing_u * wing_spn2
    wing_rel_pos = wing_z / span2

    wing_incidence_ref = (body_len + tail_len, 0, fuselage_width / 2)
    wing_incidence_ang = 2.0
    wing_forward_marching = 0.2
    wing_downward_marching = 0.8

    wing_lift_dist = EllipticLiftDist(payload, span, rho, velocity)

    wing_cl3 = np.copy(list(map(wing_lift_dist.cl_at, wing_rel_pos, wing_chord)))
    wing_cl2 = np.array([1.15 * _x for _x in wing_cl3])

    wing_lift_fig = plt.figure()
    wing_vc_ax = wing_lift_fig.add_subplot(111)
    fig_sp_n = 1000
    fig_x_sp = np.linspace(0, span2, fig_sp_n)
    fig_u_sp = np.linspace(0, 1, fig_sp_n)
    fig_vc_sp = [wing_lift_dist.velocity_circulation_at(_u) for _u in fig_u_sp]
    wing_vc_ax.plot(fig_x_sp, fig_vc_sp)
    wing_vc_ax.set_xlabel('Span-wise position')
    wing_vc_ax.set_ylabel('Velocity Circulation/({})'.format(r'$m^2 \cdot s^{-1}$'))
    wing_cl_ax = wing_vc_ax.twinx()
    wing_cl_ax.plot(wing_z, wing_cl2, '.-')
    for i in range(wing_n):
        wing_cl_ax.text(wing_z[i], wing_cl2[i], '{:.3f}'.format(wing_cl2[i]))
    wing_cl_ax.set_ylabel('Lift coefficient')
    wing_lift_fig.tight_layout()
    wing_lift_fig.show()

    wing_detail = [WingProfileParam('NACA0018', wing_chord[-13], 0, cl=wing_cl2[-13]),
                   WingProfileParam('NACA0017', wing_chord[-12], 0, cl=wing_cl2[-12]),
                   WingProfileParam('NACA0016', wing_chord[-11], 0, cl=wing_cl2[-11]),
                   WingProfileParam('NACA0015', wing_chord[-10], 0, cl=wing_cl2[-10]),
                   WingProfileParam('SC(2)-0614', wing_chord[-9], -2.523, cl=wing_cl2[-9]),
                   WingProfileParam('SC(2)-0612', wing_chord[-8], -2.180, cl=wing_cl2[-8]),
                   WingProfileParam('SC(2)-0712', wing_chord[-7], -2.637, cl=wing_cl2[-7]),
                   WingProfileParam('SC(2)-0712', wing_chord[-6], -2.439, cl=wing_cl2[-6]),
                   WingProfileParam('SC(2)-0712', wing_chord[-5], -2.303, cl=wing_cl2[-5]),
                   WingProfileParam('SC(2)-0712', wing_chord[-4], -2.418, cl=wing_cl2[-4]),
                   WingProfileParam('SC(2)-0612', wing_chord[-3], -2.306, cl=wing_cl2[-3]),
                   WingProfileParam('SC(2)-0412', wing_chord[-2], -2.071, cl=wing_cl2[-2]),
                   WingProfileParam('SC(2)-0012', wing_chord[-1], -2.000, cl=wing_cl2[-1])]

    wing_profile = construct_hwb_wing_profiles([p.airfoil for p in wing_detail],
                                               wing_chord, wing_z, wing_swp,
                                               [p.twist_ang for p in wing_detail],
                                               wing_dihedral,
                                               init_origin=wing_ref_origin,
                                               incidence=[wing_incidence_ref, wing_incidence_ang],
                                               forward=wing_forward_marching,
                                               downward=wing_downward_marching)

    '''Vertical-Stabilizer'''
    vs_root_chord = 5.5
    vs_tip_chord = 3.5
    vs_spn2 = 4.2
    vs_leading_swp = 45
    vs_planform = VSPlanform(vs_root_chord, vs_tip_chord, vs_spn2, vs_leading_swp)
    vs_u = np.array([0.00, 27.50, 55.00, 85.00, 100.00]) / 100
    vs_n = len(vs_u)
    vs_z = np.array([vs_planform.z(_u) for _u in vs_u])

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

    hs_foil = ['NACA0008'] * hs_n
    hs_cl = [hs_planform.chord_len(hs_u[i]) for i in range(hs_n)]
    hs_swp = [math.degrees(math.atan2(hs_planform.x_front(hs_u[i]), hs_z[i])) for i in range(hs_n)]

    hs_forward_marching_ratio = 0
    hs_delta_x = vs_origin[0] + vs_planform.x_front(1) - hs_forward_marching_ratio * vs_tip_chord
    hs_delta_y = vs_origin[1] + vs_spn2
    hs_origin = (hs_delta_x, hs_delta_y, 0)
    hs_profile = construct_hs_profiles(hs_foil, hs_cl, hs_z, hs_swp, origin=hs_origin)

    '''Text report'''
    print('\nLength of fuselage components:')
    print('Nose: {:.3f}\nBody: {:.3f}\nTail: {:.3f}'.format(nose_len, body_len, tail_len))
    print('\nWing:')
    print('Root chord: {:.3f}'.format(wing_root_len))
    print('Half span: {:.3f}'.format(wing_spn2))
    print('Area: {:.2f}'.format(wing_planform.area))
    print('MAC: {:.3f}'.format(wing_planform.mean_aerodynamic_chord))
    print('Dihedral: {}'.format(wing_dihedral))
    print('\nVertical-Stabilizer:')
    print('Area: {:.2f}'.format(vs_planform.area / 2))
    print('MAC: {:.3f}'.format(vs_planform.mean_aerodynamic_chord))
    print('\nHorizontal-Stabilizer:')
    print('Area: {:.2f}'.format(hs_planform.area))
    print('MAC: {:.3f}'.format(hs_planform.mean_aerodynamic_chord))

    '''Graphic representation'''
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    wing_planform.pic(ax1, u=wing_u)
    # ax2 = fig.add_subplot(111)
    # vs_planform.pic(ax2, u=vs_u, direction='vertical')
    # ax3 = fig.add_subplot(212)
    # ax3.plot(wing_u, wing_cl2, label='Cl')
    # ax3.legend()
    # ax4 = fig.add_subplot(111)
    # hs_planform.pic(ax4, u=hs_u)
    fig.tight_layout()
    plt.show()

    '''CAD model'''
    model = IGES_Model()
    # model.add(IGES_Pnt(body_front_up))
    # model.add(IGES_Pnt(body_front_down))
    # model.add(IGES_Pnt(body_front_mid))
    # model.add(IGES_Line(body_front_up, body_back_up))
    # model.add(IGES_Line(body_front_down, body_back_down))
    # model.add(IGES_Line(body_front_mid, body_back_mid))
    # model.add(body_front_crv.to_iges())
    # model.add(body_back_crv.to_iges())
    # model.add(tail_back_crv.to_iges())
    # model.add(IGES_Line(tail_back_center, tail_back_up))
    # model.add(IGES_Line(tail_back_center, tail_back_down))
    # model.add(IGES_Line(tail_back_center, tail_back_mid))
    # model.add(nose_front_crv.to_iges())
    # model.add(IGES_Line(nose_front_center, nose_front_up))
    # model.add(IGES_Line(nose_front_center, nose_front_down))
    # model.add(IGES_Line(nose_front_center, nose_front_mid))
    for _c in wing_profile:
        model.add(_c.to_iges())
    # for _c in vs_profile:
    #     model.add(_c.to_iges())
    # for _c in hs_profile:
    #     model.add(_c.to_iges())
    model.save('HWB.igs')


if __name__ == '__main__':
    construct_hwb_frame()
