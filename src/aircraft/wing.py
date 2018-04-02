from copy import deepcopy
import numpy as np
from grid import LinearTFI2D, LinearTFI3D
from grid import Plot3D
from grid.spacing import hyperbolic_tangent, single_exponential, double_exponential, uniform
from cad.iges import Model, Entity116, Entity110
from misc import pnt_pan
from cad.nurbs import Crv, Line, Spline, ConicArc
from cad.nurbs import Surf, Skinned, RuledSurf


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
