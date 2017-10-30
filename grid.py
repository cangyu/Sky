import numpy as np
import math
from scipy.linalg import norm
from copy import deepcopy
import unittest

"""
Implementation of the Plot3D standard.

Note:
All the units are SI by default.

TODO:
Import testing cases for multi-block Plot3D
"""

grid_debug = False


def report_process(msg):
    if grid_debug:
        print(msg)


class Plot3DBlock(object):
    def __init__(self, pts):
        """
        Single Block structured grid.
        IBLANK Values:
             0 : Outside of the computational area.
             1 : Normal points.
             2 : Wall boundary points.
            -n : Adjacent to the n-th block.
        :param pts: All the coordinates, the index traverse through(I, J, K) from left to right, each coordinate is consist of (X, Y, Z, IBLANK).
        """

        self.data = np.copy(pts)

    @property
    def dim(self):
        """
        Dimension of the grid.
        """

        return self.data.shape[:3]

    def __repr__(self):
        ii, jj, kk = self.dim
        return "{} x {}".format(ii, jj) if kk == 1 else "{} x {} x {}".format(ii, jj, kk)

    @property
    def pnt_num(self):
        """
        :return Number of points.
        :rtype int
        """

        ii, jj, kk = self.dim
        return ii * jj * kk

    @property
    def cell_num(self):
        """
        :return Number of cells.
        :rtype int
        """

        ii, jj, kk = self.dim
        return (ii - 1) * (jj - 1) if kk == 1 else (ii - 1) * (jj - 1) * (kk - 1)

    @property
    def face_num(self):
        """
        :return Number of faces.
        :rtype int
        """

        ii, jj, kk = self.dim
        return 3 * ii * jj * kk - (ii * jj + jj * kk + kk * ii)

    @property
    def all_pnt(self):
        t = 0
        ii, jj, kk = self.dim
        ret = np.empty((self.pnt_num, 3), float)
        for k in range(kk):
            for j in range(jj):
                for i in range(ii):
                    ret[t] = self.data[i][j][k][:3]
                    t += 1

        return ret

    def write(self, f_out, with_iblank):
        """
        Output the grid into a stream.
        :param f_out: The output stream.
        :param with_iblank: Indicate if the IBLANK info is included.
        :type with_iblank: bool
        :return: None
        """

        ii, jj, kk = self.dim
        for d in range(3):
            for k in range(kk):
                for j in range(jj):
                    for i in range(ii):
                        f_out.write("{}{}".format('\n' if i == 0 else ' ', self.data[i][j][k][d]))

        if with_iblank:
            for k in range(kk):
                for j in range(jj):
                    for i in range(ii):
                        f_out.write("{}{}".format('\n' if i == 0 else ' ', int(self.data[i][j][k][-1])))

    def set_iblank(self, i, j, k, t):
        """
        Set the IBLANK value on certain point.
        :param i: Index in X direction.
        :type i: int
        :param j: Index in Y direction.
        :type j: int
        :param k: Index in Z direction.
        :type k: int
        :param t: Target IBLANK Value.
        :type t: int
        :return: None.
        """

        self.data[i][j][k][-1] = t

    def set_area_iblank(self, rg0, rg1, rg2, t):
        """
        设置区域内网格点的IBLANK信息
        :param rg0: X方向范围
        :param rg1: Y方向范围
        :param rg2: Z方向范围
        :param t: IBLANK Value
        :type t: int
        :return: None
        """

        for i in rg0:
            for j in rg1:
                for k in rg2:
                    self.set_iblank(i, j, k, t)

    def set_boundary_iblank(self, t):
        """
        设置边界上网格点的IBLANK信息
        :param t: IBLANK Value
        :type t: int
        :return: None
        """

        ii, jj, kk = self.dim
        if kk == 1:
            for i in range(ii):
                for j in range(jj):
                    if i in (0, ii - 1) or j in (0, jj - 1):
                        self.set_iblank(i, j, 0, t)
        else:
            for i in range(ii):
                for j in range(jj):
                    for k in range(kk):
                        if i in (0, ii - 1) or j in (0, jj - 1) or k in (0, kk - 1):
                            self.set_iblank(i, j, k, t)

    @classmethod
    def construct_from_array(cls, pts):
        """
        Construct the Plot3D Block from grid array.
        :param pts: Input grid points, maybe 2 or 3 Dimensional.
                    The index traverse (I, J)/(I, J, K) from left to right,
                    Each element is consist of (X, Y, Z)/(X, Y).
        :return: Single-Block grid in Plot3D notation.
        :rtype: Plot3DBlock
        """

        if len(pts.shape) == 3:
            ni, nj, nd = pts.shape
            nk = 1
        elif len(pts.shape) == 4:
            ni, nj, nk, nd = pts.shape
        else:
            raise AssertionError("Invalid input grid array.")

        p3d = np.ones((ni, nj, nk, 4))  # Denote all the points as Normal by default.
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    for d in range(nd):
                        p3d[i][j][k][d] = pts[i][j][k][d]

        ret = cls(p3d)
        ret.set_boundary_iblank(2)  # Close the boundary by default.
        return ret


class Plot3D(object):
    def __init__(self):
        self.blk_list = []

    def __repr__(self):
        ret = "Multi-Block structured grid in Plot3D format with {} block(s)\n".format(self.size)
        for i in range(self.size):
            ret += "{}: {}\n".format(i, repr(self.blk_list[i]))

        return ret

    @property
    def size(self):
        return len(self.blk_list)

    def add(self, blk):
        """
        Append new block.
        :param blk: single block.
        :type blk: Plot3DBlock
        :return: None
        """

        self.blk_list.append(blk)

    def clear(self):
        self.blk_list.clear()

    def save(self, fn, with_iblank=False):
        """
        Output grid into file.
        :param fn: Filename.
        :type fn: str
        :param with_iblank: Indicate if the IBLANK info is included.
        :type with_iblank: bool
        :return: None
        """

        f_out = open(fn, 'w')
        f_out.write("{}".format(self.size))
        report_process("Writing grid: \'{}\' with {} block(s) ...".format(fn, self.size))

        for blk in self.blk_list:
            f_out.write("\n{} {} {}".format(blk.I, blk.J, blk.K))

        for k, blk in enumerate(self.blk_list):
            report_process("Writing block \'{}\' with dimension \'{}\' ...".format(k, repr(blk)))
            blk.write(f_out, with_iblank)

        f_out.close()
        report_process("Grid output done!")


class Plot3DTester(unittest.TestCase):
    @staticmethod
    def test_single():
        # x_min, x_max, y_min, y_max, z_min, z_max, nu, nv, nw
        rect_param = [(0, 100, 0, 60, 0, 40, 61, 16, 21),
                      (0, 100, 0, 60, 0, 40, 61, 16, 1)]
        # r_min, r_max, theta_min, theta_max, h_min, h_max, nu, nv, nw
        sect_param = [(50, 100, 60, 320, 0, 30, 61, 16, 21),
                      (50, 100, 60, 320, 0, 30, 61, 16, 1)]
        ans = []

        for p in rect_param:
            x_min, x_max, y_min, y_max, z_min, z_max, nu, nv, nw = p
            pts = np.zeros((nu, nv, nw, 3))
            u_list = np.linspace(x_min, x_max, nu)
            v_list = np.linspace(y_min, y_max, nv)
            w_list = np.linspace(z_min, z_max, nw)
            for i in range(nu):
                for j in range(nv):
                    for k in range(nw):
                        pts[i][j][k][0] = u_list[i]
                        pts[i][j][k][1] = v_list[j]
                        pts[i][j][k][2] = w_list[k]
            ans.append(pts)

        for p in sect_param:
            r_min, r_max, theta_min, theta_max, h_min, h_max, nu, nv, nw = p
            pts = np.zeros((nu, nv, nw, 3))
            u_list = np.linspace(r_min, r_max, nu)
            v_list = np.linspace(theta_min, theta_max, nv)
            w_list = np.linspace(h_min, h_max, nw)
            for i in range(nu):
                for j in range(nv):
                    for k in range(nw):
                        ct = math.radians(v_list[j])
                        pts[i][j][k] = np.array([u_list[i] * math.cos(ct), u_list[i] * math.sin(ct), w_list[k]])
            ans.append(pts)

        grid = Plot3D()
        for t in range(len(ans)):
            grid.clear()
            blk = Plot3DBlock.construct_from_array(ans[t])
            grid.add(blk)
            print(grid)
            fn = "SingleBlk-{}.xyz".format(t)
            grid.save(fn)

    def test_multi(self):
        pass


def equal_check(*args):
    if len(args) < 2:
        return True

    prev = args[0]
    for k, crd in enumerate(args):
        if k == 0:
            continue

        if not norm(prev - crd) < 1e-8:
            return False
        else:
            prev = crd

    return True


class TFI(object):
    def __init__(self):
        self.grid = None


class LinearTFI2D(TFI):
    def __init__(self, c1, c2, c3, c4):
        """
        2维无限插值(Linear)
        :param c1: 平行于x轴方向的第1条曲线，调用得到3维坐标点
        :param c2: 平行于y轴方向的第1条曲线，调用得到3维坐标点
        :param c3: 平行于x轴方向的第2条曲线，调用得到3维坐标点
        :param c4: 平行于y轴方向的第2条曲线，调用得到3维坐标点
        """

        super(LinearTFI2D, self).__init__()

        '''Defensive Check'''
        if not equal_check(c1(0), c2(0)):
            raise AssertionError("C1 start and C2 start not meet.")
        if not equal_check(c1(1), c4(0)):
            raise AssertionError("C1 end and C4 start not meet.")
        if not equal_check(c4(1), c3(1)):
            raise AssertionError("C4 end and C3 end not meet.")
        if not equal_check(c2(1), c3(0)):
            raise AssertionError("C2 end and C3 start not meet.")

        '''Copy back in case the parameters get changed outside'''
        self.c1 = deepcopy(c1)
        self.c2 = deepcopy(c2)
        self.c3 = deepcopy(c3)
        self.c4 = deepcopy(c4)

        '''Compute intersections in advance'''
        self.P12 = self.c1(0)  # c1与c2交点
        self.P14 = self.c1(1)  # c1与c4交点
        self.P23 = self.c3(0)  # c2与c3交点
        self.P34 = self.c3(1)  # c3与c4交点

        self.U = lambda u, v: (1 - u) * self.c2(v) + u * self.c4(v)
        self.V = lambda u, v: (1 - v) * self.c1(u) + v * self.c3(u)
        self.UV = lambda u, v: (1 - u) * (1 - v) * self.P12 + u * v * self.P34 + (1 - u) * v * self.P23 + (1 - v) * u * self.P14

    def __call__(self, u, v):
        """
        曲面在(u,v)处的坐标
        :param u: U方向参数
        :type u: float
        :param v: V方向参数
        :type v: float
        :return: Coordinate at (u, v)
        """

        return self.U(u, v) + self.V(u, v) - self.UV(u, v)

    def calc_grid(self, pu, pv):
        """
        根据网格点的参数分布计算对应的坐标
        :param pu: 所有网格点的U方向参数值
        :param pv: 所有网格点的V方向参数值
        :return: None
        """

        '''Dimension on each axis'''
        d1 = len(pu)
        d2 = len(pv)

        '''Complement of each parameter distribution'''
        puc = np.ones(d1) - pu
        pvc = np.ones(d2) - pv

        '''Boundary coordinates'''
        pc1 = np.empty((d1, 3), float)
        pc3 = np.empty((d1, 3), float)
        pc2 = np.empty((d2, 3), float)
        pc4 = np.empty((d2, 3), float)

        '''Utilize pre-computed results'''
        pc1[0] = self.P12
        pc1[-1] = self.P14
        pc3[0] = self.P23
        pc3[-1] = self.P34
        pc2[0] = self.P12
        pc2[-1] = self.P23
        pc4[0] = self.P14
        pc4[-1] = self.P34

        '''Compute boundary coordinates in advance'''
        for i in range(1, d1 - 1):
            u = pu[i]
            pc1[i] = self.c1(u)
            pc3[i] = self.c3(u)

        for j in range(1, d2 - 1):
            v = pv[j]
            pc2[j] = self.c2(v)
            pc4[j] = self.c4(v)

        '''Transfinite interpolation'''
        self.grid = np.empty((d1, d2, 3), float)
        for i in range(d1):
            for j in range(d2):
                self.grid[i][j] = puc[i] * pc2[j] + pu[i] * pc4[j] + pvc[j] * pc1[i] + pv[j] * pc3[i]
                self.grid[i][j] -= puc[i] * pvc[j] * self.P12 + pu[i] * pv[j] * self.P34 + puc[i] * pv[j] * self.P23 + pvc[j] * pu[i] * self.P14


class LinearTFI3D(TFI):
    def __init__(self, s1, s2, s3, s4, s5, s6):
        """
        3维无限插值(Linear),右手系
        :param s1: 垂直于x轴的第1个曲面(Left)
        :param s2: 垂直于x轴的第2个曲面(Right)
        :param s3: 垂直于y轴的第1个曲面(Front)
        :param s4: 垂直于y轴的第2个曲面(Back)
        :param s5: 垂直于z轴的第1个曲面(Bottom)
        :param s6: 垂直于z轴的第2个曲面(Top)
        """

        super(LinearTFI3D, self).__init__()

        '''Defensive check'''
        assert equal_check(s1(0, 0), s3(0, 0), s5(0, 0))
        assert equal_check(s2(0, 0), s3(0, 1), s5(1, 0))
        assert equal_check(s2(1, 0), s4(0, 1), s5(1, 1))
        assert equal_check(s1(1, 0), s4(0, 0), s5(0, 1))
        assert equal_check(s1(0, 1), s3(1, 0), s6(0, 0))
        assert equal_check(s2(0, 1), s3(1, 1), s6(1, 0))
        assert equal_check(s2(1, 1), s4(1, 1), s6(1, 1))
        assert equal_check(s1(1, 1), s4(1, 0), s6(0, 1))

        self.s1 = deepcopy(s1)
        self.s2 = deepcopy(s2)
        self.s3 = deepcopy(s3)
        self.s4 = deepcopy(s4)
        self.s5 = deepcopy(s5)
        self.s6 = deepcopy(s6)

        self.c35 = lambda u: self.s5(u, 0)  # intersection of s5, s3
        self.c15 = lambda v: self.s5(0, v)  # intersection of s5, s1
        self.c45 = lambda u: self.s5(u, 1)  # intersection of s5, s4
        self.c25 = lambda v: self.s5(1, v)  # intersection of s5, s2
        self.c36 = lambda u: self.s6(u, 0)  # intersection of s6, s3
        self.c16 = lambda v: self.s6(0, v)  # intersection of s6, s1
        self.c46 = lambda u: self.s6(u, 1)  # intersection of s6, s4
        self.c26 = lambda v: self.s6(1, v)  # intersection of s6, s2
        self.c13 = lambda w: self.s3(w, 0)  # intersection of s1, s3
        self.c23 = lambda w: self.s2(0, w)  # intersection of s3, s2
        self.c24 = lambda w: self.s4(w, 1)  # intersection of s2, s4
        self.c14 = lambda w: self.s1(1, w)  # intersection of s4, s1

        self.p135 = self.c13(0)  # intersection of s1, s3, s5
        self.p235 = self.c23(0)  # intersection of s2, s3, s5
        self.p245 = self.c24(0)  # intersection of s2, s4, s5
        self.p145 = self.c14(0)  # intersection of s1, s4, s5
        self.p136 = self.c13(1)  # intersection of s1, s3, s6
        self.p236 = self.c23(1)  # intersection of s2, s3, s6
        self.p246 = self.c24(1)  # intersection of s2, s4, s6
        self.p146 = self.c14(1)  # intersection of s1, s4, s6

        self.U = lambda u, v, w: (1 - u) * self.s1(v, w) + u * self.s2(v, w)
        self.V = lambda u, v, w: (1 - v) * self.s3(w, u) + v * self.s4(w, u)
        self.W = lambda u, v, w: (1 - w) * self.s5(u, v) + w * self.s6(u, v)
        self.UV = lambda u, v, w: (1 - u) * (1 - v) * self.c13(w) + (1 - u) * v * self.c14(w) + u * (1 - v) * self.c23(w) + u * v * self.c24(w)
        self.VW = lambda u, v, w: (1 - v) * (1 - w) * self.c35(u) + (1 - v) * w * self.c36(u) + v * (1 - w) * self.c45(u) + v * w * self.c46(u)
        self.WU = lambda u, v, w: (1 - w) * (1 - u) * self.c15(v) + (1 - w) * u * self.c25(v) + w * (1 - u) * self.c16(v) + w * u * self.c26(v)
        self.UVW = lambda u, v, w: (1 - u) * (1 - v) * (1 - w) * self.p135 + (1 - u) * (1 - v) * w * self.p136 + (1 - u) * v * (1 - w) * self.p145 + (1 - u) * v * w * self.p146 + u * (1 - v) * (1 - w) * self.p235 + u * (1 - v) * w * self.p236 + u * v * (1 - w) * self.p245 + u * v * w * self.p246

    @classmethod
    def from_edges(cls, *args):
        """
        Construct 3D TFI from its boundary edges.
        :param args: The set of boundary edges
        :return: 3D TFI
        :rtype: LinearTFI3D
        """

        assert len(args) == 12

        s1 = LinearTFI2D(args[3], args[8], args[7], args[11])
        s2 = LinearTFI2D(args[1], args[9], args[5], args[10])
        s3 = LinearTFI2D(args[8], args[0], args[9], args[4])
        s4 = LinearTFI2D(args[11], args[2], args[10], args[6])
        s5 = LinearTFI2D(args[0], args[3], args[2], args[1])
        s6 = LinearTFI2D(args[4], args[7], args[6], args[5])

        return cls(s1, s2, s3, s4, s5, s6)

    def __call__(self, u, v, w):
        """
        在(u,v,w)处的坐标
        """

        return self.U(u, v, w) + self.V(u, v, w) + self.W(u, v, w) - (self.UV(u, v, w) + self.VW(u, v, w) + self.WU(u, v, w)) + self.UVW(u, v, w)

    def calc_grid(self, pu, pv, pw):
        """
        根据网格点的参数分布计算对应的坐标
        :param pu: 所有网格点的U方向参数值
        :param pv: 所有网格点的V方向参数值
        :param pw: 所有网格点的W方向参数值
        :return: None
        """

        '''Dimension of each axis'''
        d1 = len(pu)
        d2 = len(pv)
        d3 = len(pw)

        '''Complement of each parameter distribution'''
        puc = np.ones(d1) - pu
        pvc = np.ones(d2) - pv
        pwc = np.ones(d3) - pw

        '''Boundary edge coordinates'''
        pc35 = np.empty((d1, 3), float)
        pc36 = np.empty((d1, 3), float)
        pc46 = np.empty((d1, 3), float)
        pc45 = np.empty((d1, 3), float)
        pc16 = np.empty((d2, 3), float)
        pc26 = np.empty((d2, 3), float)
        pc25 = np.empty((d2, 3), float)
        pc15 = np.empty((d2, 3), float)
        pc13 = np.empty((d3, 3), float)
        pc23 = np.empty((d3, 3), float)
        pc24 = np.empty((d3, 3), float)
        pc14 = np.empty((d3, 3), float)

        '''Boundary surface coordinates'''
        ps1 = np.empty((d2, d3, 3), float)
        ps2 = np.empty((d2, d3, 3), float)
        ps3 = np.empty((d3, d1, 3), float)
        ps4 = np.empty((d3, d1, 3), float)
        ps5 = np.empty((d1, d2, 3), float)
        ps6 = np.empty((d1, d2, 3), float)

        '''Utilize pre-computed corner coordinates'''
        pc35[0] = self.p135
        pc35[-1] = self.p235
        pc15[0] = self.p135
        pc15[-1] = self.p145
        pc45[0] = self.p145
        pc45[-1] = self.p245
        pc25[0] = self.p235
        pc25[-1] = self.p245
        pc36[0] = self.p136
        pc36[-1] = self.p236
        pc16[0] = self.p136
        pc16[-1] = self.p146
        pc46[0] = self.p146
        pc46[-1] = self.p246
        pc26[0] = self.p236
        pc26[-1] = self.p246
        pc13[0] = self.p135
        pc13[-1] = self.p136
        pc23[0] = self.p235
        pc23[-1] = self.p236
        pc24[0] = self.p245
        pc24[-1] = self.p246
        pc14[0] = self.p145
        pc14[-1] = self.p146

        '''Compute coordinates on each edge in advance'''
        for i in range(1, int(d1 - 1)):
            u = pu[i]
            pc35[i] = self.c35(u)
            pc36[i] = self.c36(u)
            pc46[i] = self.c46(u)
            pc45[i] = self.c45(u)

        for j in range(1, int(d2 - 1)):
            v = pv[j]
            pc16[j] = self.c16(v)
            pc26[j] = self.c26(v)
            pc25[j] = self.c25(v)
            pc15[j] = self.c15(v)

        for k in range(1, int(d3 - 1)):
            w = pw[k]
            pc13[k] = self.c13(w)
            pc23[k] = self.c23(w)
            pc24[k] = self.c24(w)
            pc14[k] = self.c14(w)

        '''Copy edge coordinates to surface coordinate array'''
        ps1[0] = pc13
        ps1[-1] = pc14
        ps2[0] = pc23
        ps2[-1] = pc24
        for j in range(1, int(d2 - 1)):
            ps1[j][0] = pc15[j]
            ps1[j][-1] = pc16[j]
            ps2[j][0] = pc25[j]
            ps2[j][-1] = pc26[j]

        ps3[0] = pc35
        ps3[-1] = pc36
        ps4[0] = pc45
        ps4[-1] = pc46
        for k in range(1, d3 - 1):
            ps3[k][0] = pc13[k]
            ps3[k][-1] = pc23[k]
            ps4[k][0] = pc14[k]
            ps4[k][-1] = pc24[k]

        ps5[0] = pc15
        ps5[-1] = pc25
        ps6[0] = pc16
        ps6[-1] = pc26
        for i in range(1, d1 - 1):
            ps5[i][0] = pc35[i]
            ps5[i][-1] = pc45[i]
            ps6[i][0] = pc36[i]
            ps6[i][-1] = pc46[i]

        '''Compute surface inner coordinates'''
        for j in range(1, int(d2 - 1)):
            v = pv[j]
            for k in range(1, int(d3 - 1)):
                w = pw[k]
                ps1[j][k] = self.s1(v, w)
                ps2[j][k] = self.s2(v, w)

        for k in range(1, int(d3 - 1)):
            w = pw[k]
            for i in range(1, int(d1 - 1)):
                u = pu[i]
                ps3[k][i] = self.s3(w, u)
                ps4[k][i] = self.s4(w, u)

        for i in range(1, int(d1 - 1)):
            u = pu[i]
            for j in range(1, int(d2 - 1)):
                v = pv[j]
                ps5[i][j] = self.s5(u, v)
                ps6[i][j] = self.s6(u, v)

        '''Compute the grid'''
        self.grid = np.empty((d1, d2, d3, 3), float)
        for i in range(d1):
            for j in range(d2):
                for k in range(d3):
                    self.grid[i][j][k] = puc[i] * ps1[j][k] + pu[i] * ps2[j][k]
                    self.grid[i][j][k] += pvc[j] * ps3[k][i] + pv[j] * ps4[k][i]
                    self.grid[i][j][k] += pwc[k] * ps5[i][j] + pw[k] * ps6[i][j]
                    self.grid[i][j][k] -= puc[i] * pvc[j] * pc13[k] + puc[i] * pv[j] * pc14[k] + pu[i] * pvc[j] * pc23[k] + pu[i] * pv[j] * pc24[k]
                    self.grid[i][j][k] -= pvc[j] * pwc[k] * pc35[i] + pvc[j] * pw[k] * pc36[i] + pv[j] * pwc[k] * pc45[i] + pv[j] * pw[k] * pc46[i]
                    self.grid[i][j][k] -= pwc[k] * puc[i] * pc15[j] + pwc[k] * pu[i] * pc25[j] + pw[k] * puc[i] * pc16[j] + pw[k] * pu[i] * pc26[j]
                    self.grid[i][j][k] += puc[i] * pvc[j] * pwc[k] * self.p135 + puc[i] * pvc[j] * pw[k] * self.p136 + puc[i] * pv[j] * pwc[k] * self.p145 + puc[i] * pv[j] * pw[k] * self.p146 + pu[i] * pvc[j] * pwc[k] * self.p235 + pu[i] * pvc[j] * pw[k] * self.p236 + pu[i] * pv[j] * pwc[k] * self.p245 + pu[i] * pv[j] * pw[k] * self.p246


if __name__ == '__main__':
    unittest.main()
