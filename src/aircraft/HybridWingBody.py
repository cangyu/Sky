import numpy as np
import math
from copy import deepcopy
from matplotlib import pyplot as plt
from src.iges import IGES_Model, IGES_Pnt, IGES_Line
from src.nurbs import ConicArc, LocalCubicInterpolatedCrv, point_inverse
from src.misc import pnt_pan
from src.grid import uniform
from src.wing import Wing
from src.aircraft.Baseline import global_origin, z_axis_positive, z_axis_negative
from src.aircraft.Baseline import WingPlanform, VerticalStablizerPlanform, construct_vertical_stabilizer_profiles


class HWBOuterWingPlanform(WingPlanform):
    def __init__(self, *args, **kwargs):

        cr, spn2, theta1, seg1, theta2, seg2 = args

        assert len(theta1) == len(theta2) == 3
        assert len(seg1) == len(seg2) == 2

        leading_seg_length = np.array([seg1[0], seg1[1], spn2 - sum(seg1)])
        trailing_seg_length = np.array([seg2[0], seg2[1], spn2 - sum(seg2)])
        leading_seg_theta = np.radians(theta1)
        trailing_seg_theta = np.radians(theta2)

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
        delta = (math.tan(trailing_seg_theta[0]) * trailing_seg_length[0], 0, trailing_seg_length[0])
        trailing_pnt[1] = pnt_pan(trailing_pnt[0], delta)
        if 'ct' in kwargs:
            trailing_pnt[3] = pnt_pan(leading_pnt[3], (kwargs['ct'], 0, 0))
            delta = (-math.tan(trailing_seg_theta[2]) * trailing_seg_length[2], 0, -trailing_seg_length[2])
            trailing_pnt[2] = pnt_pan(trailing_pnt[3], delta)
        else:
            for i in range(2, 4):
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


def construct_hwb_frame():
    model = IGES_Model()

    fuselage_len = 28.0
    span = 42.0
    span2 = span / 2

    fuselage_height = 4
    fuselage_width = fuselage_height * 1.2

    nose_len = 4.8
    tail_len = 7.5
    body_len = fuselage_len - (tail_len + nose_len)
    print('Length:\nNose: {:.3f}, Body: {:.3f}, Tail: {:.3f}'.format(nose_len, body_len, tail_len))

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

    # model.add(IGES_Pnt(origin))
    # model.add(IGES_Pnt(body_front_up))
    # model.add(IGES_Pnt(body_front_down))
    # model.add(IGES_Pnt(body_front_mid))
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

    fusion_width = 1.6
    outer_wing_spn2 = span2 - fuselage_width / 2 - fusion_width
    outer_wing_root_len = 0.9 * body_len + 0.1 * nose_len
    planform = HWBOuterWingPlanform(outer_wing_root_len, outer_wing_spn2,
                                    (65, 52, 28), (0.6, 1.6),
                                    (-64, -45, 15), (1.3, 3.5), ct=1.2)
    # print(planform.area)
    # print(planform.mean_aerodynamic_chord)

    u_pos = np.array([0.00, 3.33, 6.11, 12.22, 19.44, 28.33, 46.11, 82.22, 100.00]) / 100
    z_pos = u_pos * outer_wing_spn2
    n = len(u_pos)

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    planform.pic(ax1, u=u_pos)

    tc = np.array([10.0] * n) / 100
    foil = ['SC(2)-0610'] * n
    chord_len = np.array([planform.chord_len(u) for u in u_pos])
    height = chord_len * tc
    thickness_factor = np.ones(n)
    sweep_back = np.array([math.degrees(math.atan2(planform.x_front(u), planform.z(u))) for u in u_pos])
    twist = np.zeros(n)
    twist_pos = np.array([0.25] * n)
    dihedral = np.array([math.degrees(math.atan2((height[i] - height[0]) / 2, z_pos[i])) for i in range(n)])
    # dihedral = np.zeros(n)
    y_ref = np.zeros(n)

    wg = Wing.from_geom_desc(foil, chord_len, thickness_factor, z_pos, sweep_back, twist, twist_pos, dihedral, y_ref)
    for k, elem in enumerate(wg.profile):
        crv = elem.crv
        crv.pan((-1, -0.2, span2 - outer_wing_spn2))
        model.add(crv.to_iges())

    vs_root_chord = 12.5
    vs_tip_chord = 5.5
    vs_spn2 = 6.5
    vs_planform = VerticalStablizerPlanform(vs_root_chord, vs_tip_chord, vs_spn2, [1.0, 3.0], [76, 54, 36])
    vs_u = np.array([0.00, 6.67, 13.33, 27.50, 55.00, 85.00, 100.00]) / 100
    vs_n = len(vs_u)
    vs_z = np.array([vs_planform.z(_u) for _u in vs_u])

    # print(vs_planform.area / 2)
    # print(vs_planform.mean_aerodynamic_chord)

    ax2 = fig.add_subplot(122)
    vs_planform.pic(ax2, u=vs_u)

    vs_foil = ['NACA0012'] * vs_n
    vs_cl = [vs_planform.chord_len(vs_u[i]) for i in range(vs_n)]
    vs_swp = [math.degrees(math.atan2(vs_planform.x_front(vs_u[i]), vs_z[i])) for i in range(vs_n)]

    vs_delta_tail_back = 0.8
    vs_delta_body_up = -0.2
    vs_delta_x = body_len + tail_len - (vs_root_chord + vs_delta_tail_back)
    vs_delta_y = fuselage_height / 2 + vs_delta_body_up
    vs_origin = (vs_delta_x, vs_delta_y, 0)
    vs_profile = construct_vertical_stabilizer_profiles(vs_foil, vs_cl, vs_z, vs_swp, origin=vs_origin)
    for crv in vs_profile:
        model.add(crv.to_iges())

    plt.show()
    model.save('HWB2.igs')


if __name__ == '__main__':
    construct_hwb_frame()
