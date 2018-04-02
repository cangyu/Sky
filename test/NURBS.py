import unittest
import math
import numpy as np
from copy import deepcopy
from cad.iges import Model
from misc import read_airfoil_pts, sqrt2, sqrt3
from cad.nurbs import Quaternion
from cad.nurbs import find_span, all_basis_val, point_to_line, line_intersection, to_homogeneous
from cad.nurbs import Crv, Line, ConicArc, Circle, GlobalInterpolatedCrv, point_inverse, Spline
from cad.nurbs import Surf, BilinearSurf, LocalCubicInterpolatedCrv, GlobalInterpolatedSurf
from cad.nurbs import ExtrudedSurf, Coons, RuledSurf, RevolvedSurf, Skinned


class TransformTestCase(unittest.TestCase):
    def test_quaternion(self):
        # axis, angle
        data = [[(1, 1, 1), 30],
                [(1, 1, 1), 45],
                [(1, 1, 1), 90],
                [(1, 1, 1), 180]]

        # w, x, y, z
        ans = [(0.9659258, 0.1494292, 0.1494292, 0.1494292),
               (0.9238795, 0.2209424, 0.2209424, 0.2209424),
               (0.7071068, 0.4082483, 0.4082483, 0.4082483),
               (0, 0.5773503, 0.5773503, 0.5773503)]

        for k, dt in enumerate(data):
            q = Quaternion.from_u_theta(dt[0], math.radians(dt[1]))
            print(q)
            for i in range(4):
                self.assertTrue(math.isclose(ans[k][i], q.comp[i], abs_tol=1e-7))


class BasicUtilityTestCase(unittest.TestCase):
    def test_find_span(self):
        # n, p, u, u_vec
        data = [[2, 2, 0, (0, 0, 0, 1, 1, 1)],
                [2, 2, 1, (0, 0, 0, 1, 1, 1)],
                [9, 3, 0.0, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.1, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.2, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.3, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.5, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.6, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.7, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 0.8, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)],
                [9, 3, 1.0, (0, 0, 0, 0, 0.1, 0.1, 0.3, 0.5, 0.5, 0.7, 1, 1, 1, 1)]]
        ans = [2, 2, 3, 5, 5, 6, 8, 8, 9, 9, 9]

        for k, dt in enumerate(data):
            cur_ans = find_span(dt[0], dt[1], dt[2], dt[3])
            self.assertEqual(cur_ans, ans[k])

    def test_all_basis_val(self):
        # u, p, u_vec
        data = [[2.5, 0, (0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5)],
                [2.5, 1, (0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5)],
                [2.5, 2, (0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5)]]
        ans = [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
               [0, 0, 0, 0.5, 0.5, 0, 0, 0, 0],
               [0, 0, 1 / 8, 3 / 4, 1 / 8, 0, 0, 0]]

        for i in range(len(data)):
            cur_ans = all_basis_val(data[i][0], data[i][1], data[i][2])
            for j in range(len(ans[i])):
                self.assertTrue(math.isclose(cur_ans[j], ans[i][j]))

    def test_line_intersection(self):
        # p1, u1, p2, u2
        data = [[(0, 5, 0), (1, 0, 0), (0, 5, 0), (0, 0, 1)],
                [(0, 0, 0), (1, 1, 0), (5, 0, 0), (1, -1, 0)],
                [(0, 1, 0), (0, 0, 1), (0, 2, 0), (1, 0, 0)]]
        ans = [(0, 5, 0),
               (2.5, 2.5, 0),
               None]

        for i in range(len(data)):
            try:
                cur_ans = line_intersection(data[i][0], data[i][1], data[i][2], data[i][3])
            except AssertionError as e:
                print("Exception caught! with msg: \'{}\'".format(e))
            else:
                for j in range(len(ans[i])):
                    self.assertTrue(math.isclose(cur_ans[j], ans[i][j]))

    def test_point_to_line(self):
        # t, p, u
        data = [[(0, 0, 0), (5, 0, 0), (1, -1, 0)],
                [(0, 2, 0), (0, 1, 0), (0, 0, 1)],
                [(3, 4, 5), (2, 2, 2), (0, 0, 0)]]

        ans = [(2.5, 2.5, 0),
               (0, 1, 0),
               None]

        for i in range(len(ans)):
            try:
                cur_ans = point_to_line(data[i][0], data[i][1], data[i][2])
            except AssertionError as e:
                print("Exception caught! with msg: \'{}\'".format(e))
            else:
                for j in range(len(ans[i])):
                    self.assertTrue(math.isclose(cur_ans[j], ans[i][j]))


class NURBSCrvTestCase(unittest.TestCase):
    def test_construction(self):
        # knot, ctrl points, weights
        u_vec = [[0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1],
                 [0, 0, 0, 0.5, 1, 1, 1],
                 [0, 0, 0, 0, 1, 1, 1, 1],
                 [0, 0, 0, 0, 1, 1, 1, 1],
                 [0, 0, 1, 1]]
        pnt = [[[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1], [1, 0]],
               [[1, 0], [1, 1], [-1, 1], [-1, 0]],
               [[1, 0], [1, 2], [-1, 2], [-1, 0]],
               [[sqrt3 / 2, 1 / 2], [sqrt3, -3], [-sqrt3, -3], [-sqrt3 / 2, 1 / 2]],
               [[10, 10], [100, 100]]]
        w = [[1, 1 / sqrt2, 1, 1 / sqrt2, 1, 1 / sqrt2, 1, 1 / sqrt2, 1],
             [1, 0.5, 0.5, 1],
             [1, 1 / 3, 1 / 3, 1],
             [1, 1 / 6, 1 / 6, 1],
             [1, 1]]
        ans = ['circle1.igs', 'circle2.igs', 'circle3.igs', 'circle4.igs', 'line.igs']

        # Just used to avoid warning
        self.assertTrue(len(u_vec) == len(pnt) == len(w) == len(ans))

        iges_model = Model()
        for k in range(len(ans)):
            iges_model.clear()
            geom = Crv(u_vec[k], list(map(lambda _p, _w: to_homogeneous(np.append(_p, [0]), _w), pnt[k], w[k])))
            iges_model.add(geom.to_iges())
            iges_model.save(ans[k])
            print(repr(geom))

    def test_call(self):
        # knot, ctrl points, weights
        u_vec = [[0, 0, 0, 1, 2, 3, 3, 3],
                 [0, 0, 0, 1, 1, 1]]
        p = [[(0, 0, 0), (1, 1, 0), (3, 2, 0), (4, 1, 0), (5, -1, 0)],
             [(1, 0, 0), (1, 1, 0), (0, 1, 0)]]
        w = [[1, 4, 1, 1, 1],
             [1, 1, 2]]

        # u, derivative
        data = [[(0, 0), (1, 0), (3, 0)],
                [(0, 1), (0, 2), (1, 1)]]
        ans = [[(0, 0, 0), (7 / 5, 6 / 5, 0), (5, -1, 0)],
               [(0, 2, 0), (-4, 0, 0), (-1, 0, 0)]]

        # Just used to avoid warning
        self.assertTrue(len(u_vec) == len(p) == len(w) == len(data) == len(ans))

        for i in range(len(data)):
            crv = Crv(u_vec[i], list(map(lambda _p, _w: to_homogeneous(_p, _w), p[i], w[i])))
            cur_ans = list(map(lambda t: crv(t[0], t[1]), data[i]))
            np.testing.assert_array_equal(cur_ans, ans[i])

    def test_length(self):
        # knot, ctrl points, weights
        u_vec = [[0, 0, 3, 3],
                 [0, 0, 0, 1, 1, 1]]
        p = [[(0, 0, 0), (1, 1, 0)],
             [(1, 0, 0), (1, 1, 0), (0, 1, 0)]]
        w = [[1, 1],
             [1, 1, 2]]

        ans = [sqrt2, math.pi / 2]

        for i in range(len(ans)):
            crv = Crv(u_vec[i], list(map(lambda _p, _w: to_homogeneous(_p, _w), p[i], w[i])))
            cur_ans = crv.length
            self.assertTrue(math.isclose(cur_ans, ans[i]))

    def test_curvature(self):
        # knot, ctrl points, weights
        u_vec = [[0, 0, 3, 3],
                 [0, 0, 0, 1, 1, 1]]
        p = [[(0, 0, 0), (1, 1, 0)],
             [(1, 0, 0), (1, 1, 0), (0, 1, 0)]]
        w = [[1, 1],
             [1, 1, 2]]

        data = [[0, 1.5, 3],
                [0, 0.5, 1]]
        ans = [[0, 0, 0],
               [1, 1, 1]]

        for i in range(len(ans)):
            crv = Crv(u_vec[i], list(map(lambda _p, _w: to_homogeneous(_p, _w), p[i], w[i])))
            cur_ans = list(map(lambda u: crv.curvature(u), data[i]))
            for j in range(len(cur_ans)):
                self.assertTrue(math.isclose(cur_ans[j], ans[i][j]))

    def test_reverse(self):
        u_vec = [0, 0, 0, 1, 3, 6, 6, 8, 8, 8]
        s_vec = [0, 0, 0, 2, 2, 5, 7, 8, 8, 8]
        p = [(1, 2, 2), (2, 4, 8), (3, 9, 27), (4, 16, 64), (5, 25, 125), (6, 36, 216)]
        q = [(6, 36, 216), (5, 25, 125), (4, 16, 64), (3, 9, 27), (2, 4, 8), (1, 2, 2)]
        w = [1, 1, 1, 1, 1, 1]

        crv = Crv(u_vec, list(map(lambda _p, _w: to_homogeneous(_p, _w), p, w)))
        crv.reverse()
        self.assertTrue(np.array_equal(crv.U, s_vec))
        self.assertTrue(np.array_equal(crv.cpt, q))

    def test_pan(self):
        # knot, ctrl points, weights
        u_vec = [0, 0, 0, 1, 1, 1]
        p = [(1, 0, 0), (1, 1, 0), (0, 1, 0)]
        w = [1, 1, 2]

        pan_dir = (3, 4, 5)
        ans = [(4, 4, 5), (4, 5, 5), (3, 5, 5)]

        crv = Crv(u_vec, list(map(lambda _p, _w: to_homogeneous(_p, _w), p, w)))
        iges_model = Model()
        iges_model.add(crv.to_iges())
        crv.pan(pan_dir)
        iges_model.add(crv.to_iges())
        iges_model.save('test_pan.igs')

        self.assertTrue(np.array_equal(crv.cpt, ans))
        self.assertTrue(np.array_equal(crv.weight, w))

    def test_rotate(self):
        # knot, ctrl points, weights
        u_vec = [0, 0, 0, 1, 1, 1]
        p = [(1, 0, 0), (1, 1, 0), (0, 1, 0)]
        w = [1, 1, 2]

        # anchor point, rotation axis, rotation angle
        data = [[(0, 0, 0), (1, 1, 1), 30],
                [(0, 0, 0), (1, 1, 1), 45],
                [(0, 0, 0), (1, 1, 1), 90],
                [(0, 0, 0), (1, 1, 1), 180]]
        ans = [[(0.9106836, 1 / 3, -0.2440169), (2 / 3, 1.2440169, 0.0893164), (-0.2440169, 0.9106836, 1 / 3)],
               [(0.8047379, 0.5058793, -0.3106172), (0.4941207, 1.3106172, 0.1952621), (-0.3106172, 0.8047379, 0.5058793)],
               [(1 / 3, 0.9106836, -0.2440169), (0.0893164, 1.2440169, 2 / 3), (-0.2440169, 1 / 3, 0.9106836)],
               [(-1 / 3, 2 / 3, 2 / 3), (1 / 3, 1 / 3, 4 / 3), (2 / 3, -1 / 3, 2 / 3)]]

        # Just used to avoid warning
        self.assertTrue(len(data) == len(ans))

        for i, rt in enumerate(data):
            iges_model = Model()
            crv = Crv(u_vec, list(map(lambda _p, _w: to_homogeneous(_p, _w), p, w)))
            # print(crv)
            iges_model.add(crv.to_iges())
            crv.rotate(rt[0], rt[1], rt[2])
            # print(crv)
            iges_model.add(crv.to_iges())
            iges_model.save('test_rotate-{}.igs'.format(i))

            np.testing.assert_array_almost_equal(crv.cpt, ans[i])

    def test_insert_knot(self):
        u_vec = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        data = [(2.5, 1),  # 插入1个不在原节点矢量中的节点
                (2, 1),  # 插入1个在原节点矢量中的节点
                (2.5, 3),  # 插入1个不在原节点矢量中的节点3次
                (2, 2),  # 插入1个在原节点矢量中的节点2次
                (2.5, 4)]  # 插入1个不在原节点矢量中的节点4次
        ans = ['test_insert-0.igs',
               'test_insert-1.igs',
               'test_insert-2.igs',
               'test_insert-3.igs',
               None]

        for i in range(len(data)):
            crv = Crv(u_vec, pw)
            try:
                crv.insert_knot(data[i][0], data[i][1])
            except ValueError as e:
                print('Ok, illegal insertion detected with following msg:\n{}'.format(e))
            else:
                iges_model = Model()
                iges_model.add(crv.to_iges())
                iges_model.save(ans[i])

        # Should yield identical results
        crv1 = Crv(u_vec, pw)
        crv2 = Crv(u_vec, pw)
        crv1.insert_knot(2)
        crv1.insert_knot(2)
        crv2.insert_knot(2, 2)

        self.assertTrue(np.array_equal(crv1.U, crv2.U))
        self.assertTrue(np.array_equal(crv1.Pw, crv2.Pw))

    def test_refine(self):
        u_vec = [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        data = [2.5, 2.5, 2.5]
        ans = [0, 0, 0, 0, 1, 2, 2.5, 2.5, 2.5, 3, 4, 5, 5, 5, 5]

        crv = Crv(u_vec, pw)
        iges_model = Model()
        iges_model.add(crv.to_iges())
        iges_model.save('test_refine_original.igs')
        crv.refine(data)
        iges_model.clear()
        iges_model.add(crv.to_iges())
        iges_model.save('test_refine_after.igs')
        self.assertTrue(np.array_equal(crv.U, ans))

    def test_elevate(self):
        # Knot , ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 1, 1, 1, 1]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        iges_model = Model()
        crv = Crv(u_vec, pw)
        iges_model.add(crv.to_iges())
        iges_model.save('test_elevate_before.igs')
        crv.elevate(2)
        iges_model.clear()
        iges_model.add(crv.to_iges())
        iges_model.save('test_elevate_after.igs')

        self.assertTrue(True)

    def test_reparameterization(self):
        # Knot, ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.3, 0.7, 1.0, 1.0, 1.0, 1.0]
        pw = [(0, 3.14, 0, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        param = [(2, 1, 3, 2),
                 (2, 0, 3, 1)]

        iges_model = Model()
        for k, dt in enumerate(param):
            iges_model.clear()
            crv = Crv(u_vec, pw)
            iges_model.add(crv.to_iges())
            crv.reparameterization(dt[0], dt[1], dt[2], dt[3])
            iges_model.add(crv.to_iges())
            iges_model.save('test_reparameterization-{}.igs'.format(k))

        self.assertTrue(True)

    def test_split(self):
        # Knot , ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 1, 1, 1, 1]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        iges_model = Model()
        crv = Crv(u_vec, pw)
        iges_model.add(crv.to_iges())
        iges_model.save('test_split_before.igs')

        iges_model.clear()
        crv_seg = Crv.split(crv, [0.2, 0.6])
        for seg in crv_seg:
            iges_model.add(seg.to_iges())
        iges_model.save('test_split_after.igs')

        self.assertTrue(True)

    def test_decompose(self):
        # Knot , ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 1, 1, 1, 1]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]

        iges_model = Model()
        crv = Crv(u_vec, pw)
        iges_model.add(crv.to_iges())
        iges_model.save('test_decompose_before.igs')

        iges_model.clear()
        crv_seg = Crv.decompose(crv)
        for seg in crv_seg:
            iges_model.add(seg.to_iges())
        iges_model.save('test_decompose_after.igs')

        self.assertTrue(True)

    def test_point_inverse(self):
        # Knot , ctrl_pnt
        u_vec = [0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 1, 1, 1, 1]
        pw = [(0, 0, 0, 1), (0, math.pi, 0, 1), (0, 0, 4, 1), (1, 0, -2, 1), (0, 1, 0, 1), (2, 0, 0, 1), (0, 0, 9, 1), (0.618, 1.414, 2.718, 1)]
        crv = Crv(u_vec, pw)

        data = [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9, 1]
        for u in data:
            self.assertEqual(u, point_inverse(crv, crv(u)))

    def test_spline(self):
        # airfoil, degree, knot method
        data = [('M6', 3, 'chord'), ('M6', 3, 'centripetal'), ('M6', 5, 'chord'), ('M6', 5, 'centripetal'),
                ('NACA0012', 3, 'chord'), ('NACA0012', 3, 'centripetal'), ('NACA0012', 5, 'chord'), ('NACA0012', 5, 'centripetal'),
                ('RAE2822', 3, 'chord'), ('RAE2822', 3, 'centripetal'), ('RAE2822', 5, 'chord'), ('RAE2822', 5, 'centripetal')]

        iges_model = Model()
        for k, dt in enumerate(data):
            spline_foil = Spline(read_airfoil_pts(dt[0]))
            iges_model.clear()
            iges_model.add(spline_foil.to_iges())
            iges_model.save('test_spline-{}_{}_{}_{}.igs'.format(k, dt[0], dt[1], dt[2]))

        self.assertTrue(True)

    def test_global_interp(self):
        # airfoil, degree, knot method
        data = [('M6', 3, 'chord'), ('M6', 3, 'centripetal'), ('M6', 5, 'chord'), ('M6', 5, 'centripetal'),
                ('NACA0012', 3, 'chord'), ('NACA0012', 3, 'centripetal'), ('NACA0012', 5, 'chord'), ('NACA0012', 5, 'centripetal'),
                ('RAE2822', 3, 'chord'), ('RAE2822', 3, 'centripetal'), ('RAE2822', 5, 'chord'), ('RAE2822', 5, 'centripetal')]

        iges_model = Model()
        for k, dt in enumerate(data):
            crv = GlobalInterpolatedCrv(read_airfoil_pts(dt[0]), dt[1], dt[2])
            iges_model.clear()
            iges_model.add(crv.to_iges())
            iges_model.save('test_global_interp_crv-{}_{}_{}_{}.igs'.format(k, dt[0], dt[1], dt[2]))

        self.assertTrue(True)

    def test_local_interp(self):
        cr = 12
        spn = 21

        '''Front'''
        fsl = np.array([1.2, 2.3, 2.45, 0])
        fsl[-1] = spn - sum(fsl[:-1])
        alpha = np.radians([50, 45, 33, 28])
        fp = np.zeros((5, 3))
        for i in range(4):
            fp[i + 1] = fp[i]
            fp[i + 1][0] += fsl[i] * math.tan(alpha[i])
            fp[i + 1][2] += fsl[i]
        ftv = np.array([[0, 0, 1],
                        [math.sin(alpha[1]), 0, math.cos(alpha[1])],
                        [math.sin(alpha[1]), 0, math.cos(alpha[1])],
                        [math.sin(alpha[3]), 0, math.cos(alpha[3])],
                        [math.sin(alpha[3]), 0, math.cos(alpha[3])]])

        '''Tail'''
        tsl = np.array([1.6, 4.0, 1.8, 0])
        tsl[-1] = spn - sum(tsl[:-1])
        beta = np.radians([-15, -48, -20, 15])
        tp = np.zeros((5, 3))
        tp[0][0] = cr
        for i in range(4):
            tp[i + 1] = tp[i]
            tp[i + 1][0] += tsl[i] * math.tan(beta[i])
            tp[i + 1][2] += tsl[i]
        ttv = np.array([[0, 0, 1],
                        [math.sin(beta[1]), 0, math.cos(beta[1])],
                        [math.sin(beta[1]), 0, math.cos(beta[1])],
                        [math.sin(beta[3]), 0, math.cos(beta[3])],
                        [math.sin(beta[3]), 0, math.cos(beta[3])]])

        '''Display'''
        iges_model = Model()
        fc = LocalCubicInterpolatedCrv(fp, ftv)
        tc = LocalCubicInterpolatedCrv(tp, ttv)
        iges_model.add(fc.to_iges())
        iges_model.add(tc.to_iges())
        iges_model.save('lci.igs')

        self.assertTrue(True)

    def test_circle(self):
        # r, center, norm_vec
        data = [[0.3, (0, 0, 0), (0, 0, 1)],
                [0.3, (0, 0, 0), (0, 1, 0)],
                [0.3, (0, 0, 0), (1, 0, 0)],
                [0.3, (0, 0, 0), (1, 1, 1)],
                [1.5, (-1, -1, -1), (0, 0, 1)],
                [1.5, (-2.7, 0, 0), (0, 1, 0)],
                [1.5, (0, -3.14, 0), (1, 0, 0)],
                [1.5, (0, 0, -0.618), (1, 1, 1)]]

        iges_model = Model()
        for k, dt in enumerate(data):
            arc = Circle.closed_circle(dt[0], dt[1], dt[2])
            iges_model.clear()
            iges_model.add(arc.to_iges())
            iges_model.save('test_circle-{}_{}.igs'.format(k, dt[0]))

        self.assertTrue(True)

    def test_arc(self):
        # center, start, theta, norm_vec
        data = [[(0, 0, 0), (0.1, 0, 0), -20, (0, 0, 1)],
                [(0, 0, 0), (0.1, 0, 0), 0, (0, 0, 1)],
                [(0, 0, 0), (0.1, 0, 0), 3, (0, 0, 1)],
                [(0, 0, 0), (0.2, 0, 0), 45, (0, 0, 1)],
                [(0, 0, 0), (0.3, 0, 0), 65.2, (0, 0, 1)],
                [(0, 0, 0), (0.4, 0, 0), 90, (0, 0, 1)],
                [(0, 0, 0), (0.5, 0, 0), 120, (0, 0, 1)],
                [(0, 0, 0), (0.6, 0, 0), 135, (0, 0, 1)],
                [(0, 0, 0), (0.7, 0, 0), 150, (0, 0, 1)],
                [(0, 0, 0), (0.8, 0, 0), 180, (0, 0, 1)],
                [(0, 0, 0), (0.9, 0, 0), 195, (0, 0, 1)],
                [(0, 0, 0), (1.0, 0, 0), 225, (0, 0, 1)],
                [(0, 0, 0), (1.1, 0, 0), 240, (0, 0, 1)],
                [(0, 0, 0), (1.2, 0, 0), 270, (0, 0, 1)],
                [(0, 0, 0), (1.3, 0, 0), 315, (0, 0, 1)],
                [(0, 0, 0), (1.4, 0, 0), 324, (0, 0, 1)],
                [(0, 0, 0), (1.5, 0, 0), 360, (0, 0, 1)],
                [(0, 0, 0), (0.1, 0, 0), 1024, (0, 0, 1)],
                [(0.5, 1.2, 1.3), (0.1, 0, 0), 3, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.2, 0, 0), 45, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.3, 0, 0), 65.2, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.4, 0, 0), 90, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.5, 0, 0), 120, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.6, 0, 0), 135, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.7, 0, 0), 150, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.8, 0, 0), 180, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (0.9, 0, 0), 195, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.0, 0, 0), 225, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.1, 0, 0), 240, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.2, 0, 0), 270, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.3, 0, 0), 315, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.4, 0, 0), 324, (1, 1, 1)],
                [(0.5, 1.2, 1.3), (1.5, 0, 0), 360, (1, 1, 1)]]

        iges_model = Model()
        for k, dt in enumerate(data):
            try:
                arc = Circle(dt[0], dt[1], dt[2], dt[3])
            except ValueError as e:
                print('Exception caught with msg: {}'.format(e))
            else:
                iges_model.add(arc.to_iges())

        iges_model.save('test_arc.igs')
        self.assertTrue(True)

    def test_arc_2pnt(self):
        # start, end, theta, norm_vec
        data = [[(0, 0, 500), (0, 0, 100), 180, (0, 1, 0)],
                [(0, 0, 50), (0, 0, 10), 180, (0, 1, 0)],
                [(0, 0, 5), (0, 0, 1), 180, (0, 1, 0)]]

        iges_model = Model()
        for k, dt in enumerate(data):
            arc = Circle.from_2pnt(dt[0], dt[1], dt[2], dt[3])
            iges_model.clear()
            iges_model.add(arc.to_iges())
            iges_model.save('test_arc_pnt-{}.igs'.format(k))
        self.assertTrue(True)

    def test_conic(self):
        # a,b of an ellipse
        data = [(10, 6), (20, 8), (50, 12), (10, 25)]
        z = 10

        for k, dt in enumerate(data):
            a, b = dt
            arc1 = ConicArc((0, -b, 0), (1, 0, 0), (a, 0, 0), (0, 1, 0), (a / sqrt2, -b / sqrt2, 0))  # 1/4 Ellipse Circle
            arc2 = ConicArc((0, b, 0), (-1, 0, 0), (0, -b, 0), (1, 0, 0), (-a, 0, 0))  # 1/2 Ellipse Circle
            arc3 = ConicArc((0, b, z), (-1, 0, 0), (0, -b, z), (1, 0, 0), (-a, 0, z))  # 1/2 Ellipse Circle pan on Z
            iges_model = Model()
            iges_model.add(arc1.to_iges())
            iges_model.add(arc2.to_iges())
            iges_model.add(arc3.to_iges())
            print(arc2)
            iges_model.save('test_conic-{}.igs'.format(k))
        self.assertTrue(True)


class NURBSSurfTestCase(unittest.TestCase):
    def test_call(self):
        pass

    def test_bilinear(self):
        l = 10.0
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        n = len(p)
        iges_model = Model()
        for k in range(n):
            s = BilinearSurf(p[k])
            iges_model.add(s.to_iges())
        iges_model.save('test_bilinear.igs')
        self.assertTrue(True)

    def test_pan(self):
        l = 10.0
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])
        d = [(0, 0, 10), (0, 0, 10)]

        self.assertTrue(len(p) == len(d))
        n = len(p)

        iges_model = Model()
        for k in range(n):
            iges_model.clear()
            s = BilinearSurf(p[k])
            iges_model.add(s.to_iges())
            s.pan(d[k])
            iges_model.add(s.to_iges())
            iges_model.save('test_pan-{}.igs'.format(k))

    def test_rotate(self):
        l = 10.0
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])
        ref = [[l, l, l], [l, l, l]]
        axis = [[0, 0, 5], [0, 0, 5]]
        ang = [45, 45]

        n = len(p)
        iges_model = Model()
        for k in range(n):
            iges_model.clear()
            s = BilinearSurf(p[k])
            iges_model.add(s.to_iges())
            s.rotate(ref[k], axis[k], ang[k])
            iges_model.add(s.to_iges())
            iges_model.save('test_rotate-{}.igs'.format(k))
        self.assertTrue(True)

    def test_reverse(self):
        l = 12.34
        p = np.array([[[[0, 0, 0], [0, l, l]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]],
                       [[l, 0, l], [l, l, 0]]]])
        axis = ['U', 'V', 'UV', 'U', 'V', 'UV']

        iges_model = Model()
        for k, bp in enumerate(p):
            iges_model.clear()
            s = BilinearSurf(bp)
            iges_model.add(s.to_iges())
            s.reverse(axis[k])
            iges_model.add(s.to_iges())
            iges_model.save('test_reverse-{}.igs'.format(k))
        self.assertTrue(True)

    def test_swap(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        for k, bpt in enumerate(p):
            s = BilinearSurf(bpt)
            ss = deepcopy(s)
            s.elevate(1, 2)
            ss.elevate(1, 2)
            s.swap()
            iges_model.clear()
            iges_model.add(s.to_iges())
            iges_model.add(ss.to_iges())
            iges_model.save('test_swap-{}.igs'.format(k))

            ni, nj = 50, 30
            u_dist, v_dist = np.meshgrid(np.linspace(0, 1, ni), np.linspace(0, 1, nj), indexing='ij')
            ps = np.zeros((ni, nj, 3))
            pss = np.zeros((ni, nj, 3))

            for i in range(ni):
                for j in range(nj):
                    ps[i][j] = s(v_dist[i][j], u_dist[i][j])
                    pss[i][j] = ss(u_dist[i][j], v_dist[i][j])

            self.assertTrue(np.allclose(ps, pss))

    def test_insert(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        print('Test insert utility.')
        for k, dt in enumerate(p):
            print('Test-{}'.format(k))
            iges_model.clear()
            s = BilinearSurf(dt)
            print('Before:{}'.format(s))
            iges_model.add(s.to_iges())
            s.insert_knot(0.5, 1, 'U')
            s.insert_knot(0.5, 1, 'V')
            print('After:{}'.format(s))
            iges_model.add(s.to_iges())
            iges_model.save('test_insert-{}.igs'.format(k))
        self.assertTrue(True)

    def test_refine(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        print('Test refine utility.')
        for k, dt in enumerate(p):
            print('Test-{}'.format(k))
            iges_model.clear()
            s = BilinearSurf(dt)
            ss = deepcopy(s)
            print('Before:{}'.format(s))
            iges_model.add(s.to_iges())
            s.refine('U', [0.2, 0.3, 0.6, 0.7])
            s.refine('V', [0.1, 0.4, 0.8, 0.9])
            print('After:{}'.format(s))
            iges_model.add(s.to_iges())
            iges_model.save('test_refine-{}.igs'.format(k))

            ni, nj = 50, 30
            u_dist = np.linspace(0, 1, ni)
            v_dist = np.linspace(0, 1, nj)
            ps = np.zeros((ni, nj, 3))
            pss = np.zeros((ni, nj, 3))

            for i in range(ni):
                for j in range(nj):
                    ps[i][j] = s(u_dist[i], v_dist[j])
                    pss[i][j] = ss(u_dist[i], v_dist[j])

            self.assertTrue(np.allclose(ps, pss))

    def test_elevate(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        print('Test elevate utility.')
        for k, dt in enumerate(p):
            print('Test-{}'.format(k))
            iges_model.clear()
            s = BilinearSurf(dt)
            ss = deepcopy(s)
            print('Before:{}'.format(s))
            iges_model.add(s.to_iges())
            s.elevate(1, 2)
            print('After:{}'.format(s))
            iges_model.add(s.to_iges())
            iges_model.save('test_elevate-{}.igs'.format(k))

            ni, nj = 50, 30
            u_dist = np.linspace(0, 1, ni)
            v_dist = np.linspace(0, 1, nj)
            ps = np.zeros((ni, nj, 3))
            pss = np.zeros((ni, nj, 3))
            for i in range(ni):
                for j in range(nj):
                    ps[i][j] = s(u_dist[i], v_dist[j])
                    pss[i][j] = ss(u_dist[i], v_dist[j])

            self.assertTrue(np.allclose(ps, pss))

    def test_split(self):
        print('Test split utility.')
        iges_model = Model()
        s = Surf(np.array([0, 0, 0, 0, 1., 1., 1., 1.]),
                 np.array([0, 0, 0, 0.5, 1, 1, 1]),
                 np.array([[[0, 0, 0, 1.], [0, 1, 1, 1], [0, 2, 3, 1], [0, 3, 2, 1]],
                           [[1, 0, 0, 1.], [1, 1, 2, 1], [1, 3, 5, 1], [1, 4, 2, 1]],
                           [[2, 0, 0, 1.], [2, 2, 2, 1], [2, 3, 6, 1], [2, 5, 7, 1]],
                           [[3, 0, 0, 1.], [3, 1, 1, 1], [3, 2, 2, 1], [3, 3, 3, 1]]]))
        iges_model.add(s.to_iges())
        print('Original:{}'.format(s))
        ss = Surf.split(s, (0.4, 0.5, 0.6), [0.2, 0.5, 0.7])
        k = 0
        for ue in ss:
            for ve in ue:
                iges_model.add(ve.to_iges())
                print('Segment {}:{}'.format(k, ve))
                k += 1
        iges_model.save('test_split.igs')
        self.assertTrue(True)

    def test_extract(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])
        ext = [('U', 0), ('U', 0.5), ('U', 1), ('V', 0), ('V', 0.5), ('V', 1),
               ('U', 0), ('U', 0.5), ('U', 1), ('V', 0), ('V', 0.5), ('V', 1)]
        self.assertTrue(len(p) == len(ext))

        iges_model = Model()
        print('Test extract utility.')
        for k, dt in enumerate(p):
            print('Case-{}'.format(k))
            s = BilinearSurf(dt)
            print('Surf:\n{}'.format(s))
            mid = s.extract(ext[k][0], ext[k][1])
            print('Curve at {}={}:\n{}'.format(ext[k][0], ext[k][1], mid))
            iges_model.clear()
            iges_model.add(s.to_iges())
            iges_model.add(mid.to_iges())
            iges_model.save('test_extract-{}_{}={}.igs'.format(k, ext[k][0], ext[k][1]))

    def test_extrude(self):
        crv = [Circle.from_2pnt((30, 45, 69), (44, 66, 88), 75, (1, -1, 1)),
               Circle.from_2pnt((0, 50, 0), (50, 0, 0), 315, (0, 0, 1))]
        direction = [(30, 30, 30), (0, 0, 90)]
        self.assertTrue(len(crv) == len(direction))

        iges_model = Model()
        for k in range(len(crv)):
            sf = ExtrudedSurf(crv[k], direction[k])
            iges_model.clear()
            iges_model.add(sf.to_iges())
            iges_model.save('test_extrude-{}.igs'.format(k))

    def test_revolved(self):
        generatrix = [Line((50, 0, 0), (50, 0, 20)),
                      Circle.from_2pnt((100, 0, 0), (100, 0, 20), 180, (1, 0, 0))]
        center = [(20, 0, 0), (50, 0, 0)]
        axis = [(0, 0, 1), (0, 0, 1)]
        theta = [60, 180]

        iges_model = Model()
        for k in range(len(generatrix)):
            iges_model.clear()
            iges_model.add(generatrix[k].to_iges())
            sf = RevolvedSurf(center[k], axis[k], theta[k], generatrix[k])
            iges_model.add(sf.to_iges())
            iges_model.save('test_revolved-{}.igs'.format(k))
        self.assertTrue(True)

    def test_ruled(self):
        # foil, z
        data = [('M6', 0), ('M6', 0.1), ('M6', 0.2), ('M6', 0.3),
                ('NACA0012', 0.4), ('NACA0012', 0.5), ('NACA0012', 0.6), ('NACA0012', 0.7),
                ('RAE2822', 0.8), ('RAE2822', 0.9), ('RAE2822', 1.0), ('RAE2822', 1.1)]
        crv = []
        for k in range(len(data)):
            foil = data[k][0]
            z = data[k][1]
            pts = read_airfoil_pts(foil)
            c = Spline(pts)
            c.pan((0, 0, z))
            crv.append(c)
        self.assertTrue(len(crv) == len(data))
        iges_model = Model()
        for i in range(1, len(crv)):
            sf = RuledSurf(crv[i - 1], crv[i])
            iges_model.add(sf.to_iges())
        iges_model.save('test_ruled.igs')

    def test_global_interp(self):
        # foil, N: num of profile, L: length per profile, p, q
        data = [('M6', 10, 0.7, 3, 3),
                ('M6', 10, 0.7, 3, 5),
                ('M6', 10, 0.7, 5, 3),
                ('M6', 10, 0.7, 5, 5),
                ('NACA0012', 10, 0.7, 3, 3),
                ('NACA0012', 10, 0.7, 3, 5),
                ('NACA0012', 10, 0.7, 5, 3),
                ('NACA0012', 10, 0.7, 5, 5),
                ('RAE2822', 10, 0.7, 3, 3),
                ('RAE2822', 10, 0.7, 3, 5),
                ('RAE2822', 10, 0.7, 5, 3),
                ('RAE2822', 10, 0.7, 5, 5)]

        iges_model = Model()
        for k in range(len(data)):
            foil = data[k][0]
            n = data[k][1]
            l = data[k][2]
            p = data[k][3]
            q = data[k][4]
            pts = read_airfoil_pts(foil)
            m, dim = pts.shape
            all_pts = np.zeros((n + 1, m, dim))
            for i in range(n + 1):
                all_pts[i] = np.copy(pts)
                for j in range(m):
                    all_pts[i][j][-1] = l * i
            wsf = GlobalInterpolatedSurf(all_pts, p, q)
            iges_model.clear()
            iges_model.add(wsf.to_iges())
            iges_model.save("test_global_interp-{}_{}_{}_{}.igs".format(k, foil, p, q))
        self.assertTrue(True)

    def test_local_interp(self):
        pass

    def test_coons(self):
        l = 20
        u0 = Spline(read_airfoil_pts('M6') * l)
        u0.pan((0, 0, 5))
        u1 = Circle.from_2pnt((25, 60, 5), (25, -60, 5), 180, (0, 0, 1))
        v0 = Line(u0.start, u1.start)
        v1 = Line(u0.end, u1.end)
        s = Coons(u0, u1, v0, v1)

        iges_model = Model()
        iges_model.add(s.to_iges())
        iges_model.add(u0.to_iges())
        iges_model.add(u1.to_iges())
        iges_model.add(v0.to_iges())
        iges_model.add(v1.to_iges())
        iges_model.save('test_coons.igs')
        self.assertTrue(True)

    def test_skinned(self):
        # foil, N: num of profile, L: length per profile, p, q
        data = [('M6', 10, 0.7, 3, 3),
                ('M6', 10, 0.7, 3, 5),
                ('M6', 10, 0.7, 5, 3),
                ('M6', 10, 0.7, 5, 5),
                ('NACA0012', 10, 0.7, 3, 3),
                ('NACA0012', 10, 0.7, 3, 5),
                ('NACA0012', 10, 0.7, 5, 3),
                ('NACA0012', 10, 0.7, 5, 5),
                ('RAE2822', 10, 0.7, 3, 3),
                ('RAE2822', 10, 0.7, 3, 5),
                ('RAE2822', 10, 0.7, 5, 3),
                ('RAE2822', 10, 0.7, 5, 5)]

        iges_model = Model()
        for k in range(len(data)):
            foil = data[k][0]
            n = data[k][1]
            l = data[k][2]
            p = data[k][3]
            q = data[k][4]
            pts = read_airfoil_pts(foil)
            m, dim = pts.shape
            all_pts = np.zeros((n + 1, m, dim))
            for i in range(n + 1):
                all_pts[i] = np.copy(pts)
                for j in range(m):
                    all_pts[i][j][-1] = l * i

            crv_list = [GlobalInterpolatedCrv(all_pts[i], p) for i in range(n + 1)]
            wsf = Skinned(crv_list, p, q)
            iges_model.clear()
            iges_model.add(wsf.to_iges())
            iges_model.save("test_Skinned-{}_{}_{}_{}.igs".format(k, foil, p, q))
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
