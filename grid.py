import unittest
import sys
import math
import numpy as np
import pyamg
from copy import deepcopy
from abc import abstractmethod
from scipy.linalg import norm
from scipy import sparse
from scipy.sparse.linalg import dsolve

"""
Implementation of the Plot3D standard.

Note:
All the units are SI by default.
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
        :param pts: Input grid points, the index traverse (I, J)/(I, J, K) from left to right,
                    point should be 3D, each element is consist of (X, Y, Z).
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

        p3d = np.zeros((ni, nj, nk, 4))
        if len(pts.shape) == 3:
            for i in range(ni):
                for j in range(nj):
                    for d in range(nd):
                        p3d[i][j][0][d] = pts[i][j][d]
                    p3d[i][j][0][-1] = 1
        else:
            for i in range(ni):
                for j in range(nj):
                    for k in range(nk):
                        for d in range(nd):
                            p3d[i][j][k][d] = pts[i][j][k][d]
                        p3d[i][j][k][-1] = 1

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
            ii, jj, kk = blk.dim
            f_out.write("\n{} {} {}".format(ii, jj, kk))

        for k, blk in enumerate(self.blk_list):
            report_process("->Writing block \'{}\' with dimension \'{}\' ...".format(k, repr(blk)))
            blk.write(f_out, with_iblank)

        f_out.close()
        report_process("Grid output done!")


class Plot3DTester(unittest.TestCase):
    def test_single(self):
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

        self.assertTrue(len(ans) == len(rect_param) + len(sect_param))

        grid = Plot3D()
        for t in range(len(ans)):
            grid.clear()
            blk = Plot3DBlock.construct_from_array(ans[t])
            grid.add(blk)
            fn = "'test_plot3d_single-{}.xyz".format(t)
            grid.save(fn)

    def test_multi(self):
        # x_min, x_max, y_min, y_max, z_min, z_max, nu, nv, nw
        rect_param = [(0, 100, 0, 60, 0, 40, 61, 16, 21),
                      (120, 200, 75, 100, 50, 80, 61, 16, 21)]
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

        self.assertTrue(len(ans) == len(rect_param))

        grid = Plot3D()
        for t in range(len(ans)):
            blk = Plot3DBlock.construct_from_array(ans[t])
            grid.add(blk)
        grid.save('test_plot3d_multi.xyz')


"""
Implementation of the Trans-finite Interpolation(TFI) technique.

Note:
Here only linear version are implemented.
"""


def equal_check(*args):
    if len(args) < 2:
        return True

    prev = args[0]
    for k in range(1, len(args)):
        crd = args[k]
        if not norm(prev - crd) < 1e-8:
            return False
        else:
            prev = crd

    return True


def share(p, a, b):
    return (1 - p) * a + p * b


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


class TFITester(unittest.TestCase):
    def test_2d_rect(self):
        # L, W, U, V
        data = [(5, 4, 11, 9),
                (8, 8, 31, 21)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([dt[0] * u, 0, 0]),  # C1
                              lambda v: np.array([0, dt[1] * v, 0]),  # C3
                              lambda u: np.array([dt[0] * u, dt[1], 0]),  # C2
                              lambda v: np.array([dt[0], dt[1] * v, 0]))  # C4
            tfi.calc_grid(np.linspace(0, 1, dt[2]),
                          np.linspace(0, 1, dt[3]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_rect-{}.xyz'.format(k))

    def test_2d_circle(self):
        # R1, R2, U, V
        data = [(1, 2, 6, 11),
                (0, 5, 16, 33)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([(1 - u) * dt[0] + u * dt[1], 0, 0]),
                              lambda v: np.array([dt[0] * math.cos(0.5 * math.pi * v), dt[0] * math.sin(0.5 * math.pi * v), 0]),
                              lambda u: np.array([0, (1 - u) * dt[0] + u * dt[1], 0]),
                              lambda v: np.array([dt[1] * math.cos(0.5 * math.pi * v), dt[1] * math.sin(0.5 * math.pi * v), 0]))

            tfi.calc_grid(np.linspace(0, 1, dt[2]),
                          np.linspace(0, 1, dt[3]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_circle-{}.xyz'.format(k))

    def test_2d_eccentric(self):
        # Delta, R1, R2, U, V
        data = [(-10, 4, 25, 16, 41),
                (-10, 0, 25, 16, 41),
                (-1, 2, 5, 16, 33)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([(dt[0] + dt[1]) * (1 - u) + dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[1] * math.cos(math.pi * v) + dt[0], dt[1] * math.sin(math.pi * v), 0]),
                              lambda u: np.array([(dt[0] - dt[1]) * (1 - u) - dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[2] * math.cos(math.pi * v), dt[2] * math.sin(math.pi * v), 0]))

            tfi.calc_grid(np.linspace(0, 1, dt[3]),
                          np.linspace(0, 1, dt[4]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_eccentric-{}.xyz'.format(k))

    def test_2d_crv_rect(self):
        # L, H1, H2, H3
        data = [(100, 40, 60, 10, 50, 25)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([u * dt[0], 4 * dt[3] * u * (1 - u), 0]),
                              lambda v: np.array([0, v * dt[1], 0]),
                              lambda u: np.array([u * dt[0], (dt[1] * (1 - u * u) + dt[2] * u * u), 0]),
                              lambda v: np.array([dt[0], v * dt[2], 0]))

            tfi.calc_grid(np.linspace(0, 1, dt[4]),
                          np.linspace(0, 1, dt[5]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_crv_rect-{}.xyz'.format(k))

    def test_3d_cuboid(self):
        # L, W, H, U, V, W
        data = [(5, 4, 3, 11, 9, 5),
                (10, 10, 10, 21, 31, 41)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI3D(lambda v, w: np.array([0, v * dt[1], w * dt[2]]),
                              lambda v, w: np.array([dt[0], v * dt[1], w * dt[2]]),
                              lambda w, u: np.array([u * dt[0], 0, w * dt[2]]),
                              lambda w, u: np.array([u * dt[0], dt[1], w * dt[2]]),
                              lambda u, v: np.array([u * dt[0], v * dt[1], 0]),
                              lambda u, v: np.array([u * dt[0], v * dt[1], dt[2]]))
            tfi.calc_grid(np.linspace(0, 1.0, dt[3]),
                          np.linspace(0, 1.0, dt[4]),
                          np.linspace(0, 1.0, dt[5]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_3d_cuboid-{}.xyz'.format(k))

    def test_3d_sect(self):
        # R_MIN, R_MAX, THETA_MIN, THETA_MAX, H_MIN, H_MAX, U, V, W
        data = [(5, 20, 60, 120, -50, 50, 31, 16, 61)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI3D(lambda v, w: np.array([dt[0] * math.cos(math.radians(share(v, dt[2], dt[3]))), dt[0] * math.sin(math.radians(share(v, dt[2], dt[3]))), share(w, dt[4], dt[5])]),
                              lambda v, w: np.array([dt[1] * math.cos(math.radians(share(v, dt[2], dt[3]))), dt[1] * math.sin(math.radians(share(v, dt[2], dt[3]))), share(w, dt[4], dt[5])]),
                              lambda w, u: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(dt[2])), share(u, dt[0], dt[1]) * math.sin(math.radians(dt[2])), share(w, dt[4], dt[5])]),
                              lambda w, u: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(dt[3])), share(u, dt[0], dt[1]) * math.sin(math.radians(dt[3])), share(w, dt[4], dt[5])]),
                              lambda u, v: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(share(v, dt[2], dt[3]))), share(u, dt[0], dt[1]) * math.sin(math.radians(share(v, dt[2], dt[3]))), dt[4]]),
                              lambda u, v: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(share(v, dt[2], dt[3]))), share(u, dt[0], dt[1]) * math.sin(math.radians(share(v, dt[2], dt[3]))), dt[5]]))
            tfi.calc_grid(np.linspace(0, 1.0, dt[6]),
                          np.linspace(0, 1.0, dt[7]),
                          np.linspace(0, 1.0, dt[8]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_3d_sect-{}.xyz'.format(k))


"""
Implementation of the grid smoothing tools using elliptic PDE operator.

Note:
1. (i,j,k) is corresponding to (x,y,z),(u,v.w),(xi, eta, zeta),(x1,x2,x3),(I,J,K)...
2. (N,M) notation is not suggested to avoid confuse on column-and-row
3. By default, delta_xi = delta_eta = delta_zeta = 1
"""


class EllipticGrid2D(object):
    def __init__(self, grid):
        """
        二维曲边结构网格
        :param grid: 二维初始网格点数组，下标依次迭代I, J. 每个元素包含(X, Y, Z) / (X, Y)
        """

        if len(grid.shape) != 3:
            raise AssertionError("Invalid input grid!")

        self.I, self.J, self.Dim = grid.shape
        self.r = np.zeros((self.I, self.J, 3))  # 不论输入数组中每个元素是(X, Y)还是(X, Y, Z)，均转换成3维存储
        for i in range(self.I):
            for j in range(self.J):
                for d in range(self.Dim):
                    self.r[i][j][d] = grid[i][j][d]

        self.alpha = np.zeros((self.I, self.J))
        self.beta = np.zeros((self.I, self.J))
        self.gamma = np.zeros((self.I, self.J))

        self.I -= 1  # 网格点在第1维度的最大下标
        self.J -= 1  # 网格点在第2维度的最大下标

        self.delta_xi = 1.0
        self.delta_eta = 1.0
        self.eps = 1e-5  # 残差限

    def set_step(self, xi, eta):
        self.delta_xi = xi
        self.delta_eta = eta

    def set_eps(self, eps):
        self.eps = eps

    def get_grid(self):
        return self.r

    @abstractmethod
    def coefficient_around(self, i, j):
        pass

    @classmethod
    def coordinate_around(cls, i, j):
        """
        (i,j)周围的9个点的坐标
        """

        return np.array([(i, j), (i + 1, j), (i, j + 1), (i - 1, j), (i, j - 1), (i + 1, j + 1), (i - 1, j + 1), (i - 1, j - 1), (i + 1, j - 1)])

    def index_around(self, k):
        """
        第k个点周围9个点的序号
        """

        return np.array([k, k + 1, k + self.I - 1, k - 1, k - self.I + 1, k + self.I, k + self.I - 2, k - self.I, k - self.I + 2])

    def is_special(self, i, j):
        """
        判断坐标是否在边界上
        """

        return True if i in (0, self.I) or j in (0, self.J) else False

    def pder(self, i, j, comp, dir1, dir2=None):
        """
        中心差分偏导数
        :param i: 网格点在第1维度的下标
        :param j: 网格点在第2维度的下标
        :param comp: 分量
        :param dir1: 偏导数方向
        :param dir2: 偏导数方向
        :return: 在(i,j)处的偏导数
        """

        if comp not in ('x', 'y'):
            raise AssertionError("Invalid component!")
        if dir1 not in ('xi', 'eta'):
            raise AssertionError("Invalid first direction!")
        if (dir2 is not None) and (dir2 not in ('xi', 'eta')):
            raise AssertionError("Invalid second direction!")

        comp = ord(comp) - ord('x')

        if dir2 is None:
            if dir1 == 'xi':
                return (self.r[i + 1][j][comp] - self.r[i - 1][j][comp]) / (2.0 * self.delta_xi)
            else:
                return (self.r[i][j + 1][comp] - self.r[i][j - 1][comp]) / (2.0 * self.delta_eta)
        else:
            if dir1 == dir2:
                if dir1 == 'xi':
                    return (self.r[i + 1][j][comp] - 2.0 * self.r[i][j][comp] + self.r[i - 1][j][comp]) / (self.delta_xi ** 2)
                else:
                    return (self.r[i][j + 1][comp] - 2.0 * self.r[i][j][comp] + self.r[i][j - 1][comp]) / (self.delta_eta ** 2)
            else:
                return (self.r[i + 1][j + 1][comp] + self.r[i - 1][j - 1][comp] - self.r[i + 1][j - 1][comp] - self.r[i - 1][j + 1][comp]) / (4.0 * self.delta_xi * self.delta_eta)

    def pvec(self, i, j, dir1, dir2=None):
        """
        中心差分偏导矢
        :param i: 网格点在第1维度的下标
        :param j: 网格点在第2维度的下标
        :param dir1: 偏导数方向
        :param dir2: 偏导数方向
        :return: 在(i,j)处的偏导数
        """

        if dir1 not in ('xi', 'eta'):
            raise AssertionError("Invalid first direction!")
        if (dir2 is not None) and (dir2 not in ('xi', 'eta')):
            raise AssertionError("Invalid second direction!")

        if dir2 is None:
            if dir1 == 'xi':
                return (self.r[i + 1][j] - self.r[i - 1][j]) / (2.0 * self.delta_xi)
            else:
                return (self.r[i][j + 1] - self.r[i][j - 1]) / (2.0 * self.delta_eta)
        else:
            if dir1 == dir2:
                if dir1 == 'xi':
                    return (self.r[i + 1][j] - 2.0 * self.r[i][j] + self.r[i - 1][j]) / (self.delta_xi ** 2)
                else:
                    return (self.r[i][j + 1] - 2.0 * self.r[i][j] + self.r[i][j - 1]) / (self.delta_eta ** 2)
            else:
                return (self.r[i + 1][j + 1] - self.r[i + 1][j - 1] - self.r[i - 1][j + 1] + self.r[i - 1][j - 1]) / (4.0 * self.delta_xi * self.delta_eta)

    def calc_alpha_beta_gamma(self):
        """
        计算所有网格点上的alpha，beta, gamma
        """

        for i in range(1, self.I):
            for j in range(1, self.J):
                xx = self.pder(i, j, 'x', 'xi')
                xe = self.pder(i, j, 'x', 'eta')
                yx = self.pder(i, j, 'y', 'xi')
                ye = self.pder(i, j, 'y', 'eta')
                self.alpha[i][j] = xe ** 2 + ye ** 2
                self.beta[i][j] = xx * xe + yx * ye
                self.gamma[i][j] = xx ** 2 + yx ** 2

    def calc_coefficient_matrix(self):
        """
        Calculate the coefficient matrix for the 9 point stencil
        """

        var_cnt = (self.I - 1) * (self.J - 1)
        b = np.zeros((2, var_cnt))
        rows = []
        cols = []
        val = []
        ck = 0  # equation counter

        for j in range(1, self.J):
            for i in range(1, self.I):
                a = self.coefficient_around(i, j)
                coord = self.coordinate_around(i, j)
                index = self.index_around(ck)
                for k in range(9):
                    ii, jj = coord[k]
                    if self.is_special(ii, jj):
                        b[0][ck] -= a[k] * self.r[ii][jj][0]
                        b[1][ck] -= a[k] * self.r[ii][jj][1]
                    else:
                        rows.append(ck)
                        cols.append(index[k])
                        val.append(a[k])

                ck += 1

        '''Construct sparse coefficient matrix'''
        scm = sparse.coo_matrix((val, (rows, cols)), shape=(var_cnt, var_cnt), dtype=float)
        scm_csr = scm.tocsr()
        return scm_csr, b

    def calc_grid(self, choice='Picard', w=1.0, use_amg=False):
        """
        SOR迭代次数较多，求解较慢
        Picard迭代基于9点格式，通过求解稀疏矩阵的方式来求解网格，求解较快
        :return: None
        """

        if choice == 'SOR':
            if w is None:
                raise ValueError('Should specify relax factor!')
            else:
                self.sor_iteration(w)
        elif choice == 'Picard':
            if use_amg is None:
                raise ValueError('Should specify whether using a AMG solver！')
            else:
                self.picard_iteration(use_amg)
        else:
            raise ValueError('Invalid solution choice!')

    @abstractmethod
    def sor_iteration(self, w):
        pass

    @abstractmethod
    def picard_iteration(self, use_amg):
        pass

    @classmethod
    def calc_diff(cls, grid1, grid2):
        """
        计算两套网格间差值
        """

        if grid1.shape != grid2.shape:
            raise AssertionError("Grid shape do not match!")

        tmp = grid1 - grid2
        ans = 0.0
        for i in range(len(tmp)):
            for j in range(len(tmp[i])):
                ans = max(ans, np.linalg.norm(tmp[i][j], np.inf))

        return ans


class Laplace2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        对网格做Laplace光顺
        :param grid: 初始网格
        """

        super(Laplace2D, self).__init__(grid)

    def sor_iteration(self, w):
        """
        超松弛迭代(SOR)计算内部网格点
        :param w: 松弛因子
        :return: None
        """

        iteration_cnt = 0
        residual = sys.float_info.max
        while residual > self.eps:
            self.calc_alpha_beta_gamma()
            new_grid = np.copy(self.r)

            for i in range(1, self.I):
                for j in range(1, self.J):
                    t1 = self.alpha[i][j] / (self.delta_xi ** 2) * (self.r[i + 1][j] + self.r[i - 1][j])
                    t2 = 2.0 * self.beta[i][j] * np.array([self.pder(i, j, 'x', 'xi', 'eta'), self.pder(i, j, 'y', 'xi', 'eta'), 0])
                    t3 = self.gamma[i][j] / (self.delta_eta ** 2) * (self.r[i][j + 1] + self.r[i][j - 1])
                    t0 = 2.0 * (self.alpha[i][j] / (self.delta_xi ** 2) + self.gamma[i][j] / (self.delta_eta ** 2))
                    t = (t1 - t2 + t3) / t0
                    new_grid[i][j] = w * t + (1.0 - w) * self.r[i][j]

            iteration_cnt += 1
            residual = EllipticGrid2D.calc_diff(self.r, new_grid)
            print("{}: {}".format(iteration_cnt, residual))
            self.r = np.copy(new_grid)

    def coefficient_around(self, i, j):
        ax2 = self.alpha[i][j] / self.delta_xi ** 2
        ge2 = self.gamma[i][j] / self.delta_eta ** 2
        bxe2 = self.beta[i][j] / (2.0 * self.delta_xi * self.delta_eta)

        return np.array([-2 * (ax2 + ge2), ax2, ge2, ax2, ge2, -bxe2, bxe2, -bxe2, bxe2], float)

    def picard_iteration(self, use_amg):
        iteration_cnt = 0
        residual = sys.float_info.max
        while residual > self.eps:
            u = np.zeros((2, (self.I - 1) * (self.J - 1)))
            old_grid = np.copy(self.r)

            self.calc_alpha_beta_gamma()
            A, b = self.calc_coefficient_matrix()

            if use_amg:
                ml = pyamg.ruge_stuben_solver(A)
                u[0] = ml.solve(b[0], tol=1e-10)
                u[1] = ml.solve(b[1], tol=1e-10)
            else:
                u[0] = dsolve.spsolve(A, b[0], use_umfpack=True)
                u[1] = dsolve.spsolve(A, b[1], use_umfpack=True)

            cnt = 0
            for j in range(1, self.J):
                for i in range(1, self.I):
                    self.r[i][j][0] = u[0][cnt]
                    self.r[i][j][1] = u[1][cnt]
                    cnt += 1

            iteration_cnt += 1
            residual = EllipticGrid2D.calc_diff(self.r, old_grid)
            print("{}: {}".format(iteration_cnt, residual))


class ThomasMiddlecoff2D(EllipticGrid2D):
    def __init__(self, grid):
        """
        基于Thomas-Middlecoff方法对网格进行光顺
        """

        super(ThomasMiddlecoff2D, self).__init__(grid)
        self.phi = np.zeros((self.I + 1, self.J + 1))
        self.psi = np.zeros((self.I + 1, self.J + 1))

    def calc_phi_psi(self):
        """
        计算所有网格点上的phi, psi
        """

        '''Boundary'''
        for j in (0, self.J):
            for i in range(1, self.I):
                xx = self.pder(i, j, 'x', 'xi')
                yx = self.pder(i, j, 'y', 'xi')
                xxx = self.pder(i, j, 'x', 'xi', 'xi')
                yxx = self.pder(i, j, 'y', 'xi', 'xi')
                self.phi[i][j] = -(xx * xxx + yx * yxx) / (xx ** 2 + yx ** 2)

        for i in (0, self.I):
            for j in range(1, self.J):
                xe = self.pder(i, j, 'x', 'eta')
                ye = self.pder(i, j, 'y', 'eta')
                xee = self.pder(i, j, 'x', 'eta', 'eta')
                yee = self.pder(i, j, 'y', 'eta', 'eta')
                self.psi[i][j] = -(xe * xee + ye * yee) / (xe ** 2 + ye ** 2)

        '''Linear Interpolate'''
        for i in range(1, self.I):
            dist = np.linspace(self.phi[i][0], self.phi[i][self.J], self.J + 1)
            for j in range(1, self.J):
                self.phi[i][j] = dist[j]

        for j in range(1, self.J):
            dist = np.linspace(self.psi[0][j], self.psi[self.I][j], self.I + 1)
            for i in range(1, self.I):
                self.psi[i][j] = dist[i]

    def sor_iteration(self, w):
        """
        超松弛迭代计算内部网格点
        :param w: 松弛因子
        :return: None
        """

        iteration_cnt = 0
        residual = sys.float_info.max
        while residual > self.eps:
            self.calc_alpha_beta_gamma()
            self.calc_phi_psi()
            new_grid = np.copy(self.r)

            for i in range(1, self.I):
                for j in range(1, self.J):
                    t1 = self.alpha[i][j] * (1.0 / self.delta_xi ** 2 + self.phi[i][j] / (2 * self.delta_xi)) * self.r[i + 1][j]
                    t2 = self.alpha[i][j] * (1.0 / self.delta_xi ** 2 - self.phi[i][j] / (2 * self.delta_xi)) * self.r[i - 1][j]
                    t3 = 2 * self.beta[i][j] * self.pvec(i, j, 'xi', 'eta')
                    t4 = self.gamma[i][j] * (1.0 / self.delta_eta ** 2 + self.psi[i][j] / (2 * self.delta_eta)) * self.r[i][j + 1]
                    t5 = self.gamma[i][j] * (1.0 / self.delta_eta ** 2 - self.psi[i][j] / (2 * self.delta_eta)) * self.r[i][j - 1]
                    t0 = 2 * (self.alpha[i][j] / self.delta_xi ** 2 + self.gamma[i][j] / self.delta_eta ** 2)
                    t = (t1 + t2 - t3 + t4 + t5) / t0
                    new_grid[i][j] = w * t + (1 - w) * self.r[i][j]

            iteration_cnt += 1
            residual = EllipticGrid2D.calc_diff(self.r, new_grid)
            print("{}: {}".format(iteration_cnt, residual))
            self.r = np.copy(new_grid)

    def coefficient_around(self, i, j):
        ax2 = self.alpha[i][j] / self.delta_xi ** 2
        ge2 = self.gamma[i][j] / self.delta_eta ** 2
        bxe2 = self.beta[i][j] / (2.0 * self.delta_xi * self.delta_eta)
        x21 = 1.0 / self.delta_xi ** 2
        e21 = 1.0 / self.delta_eta ** 2
        phx2 = self.phi[i][j] / (2.0 * self.delta_xi)
        pse2 = self.psi[i][j] / (2.0 * self.delta_eta)

        return np.array([-2.0 * (ax2 + ge2),
                         self.alpha[i][j] * (x21 + phx2),
                         self.gamma[i][j] * (e21 + pse2),
                         self.alpha[i][j] * (x21 - phx2),
                         self.gamma[i][j] * (e21 - pse2),
                         -bxe2, bxe2, -bxe2, bxe2])

    def picard_iteration(self, use_amg):
        """
        9点格式矩阵迭代计算内部网格点
        """

        iteration_cnt = 0
        residual = sys.float_info.max
        while residual > self.eps:
            u = np.zeros((2, (self.I - 1) * (self.J - 1)))
            old_grid = np.copy(self.r)

            self.calc_alpha_beta_gamma()
            self.calc_phi_psi()
            A, b = self.calc_coefficient_matrix()

            if use_amg:
                ml = pyamg.ruge_stuben_solver(A)
                u[0] = ml.solve(b[0], tol=1e-10)
                u[1] = ml.solve(b[1], tol=1e-10)
            else:
                u[0] = dsolve.spsolve(A, b[0], use_umfpack=True)
                u[1] = dsolve.spsolve(A, b[1], use_umfpack=True)

            cnt = 0
            for j in range(1, self.J):
                for i in range(1, self.I):
                    self.r[i][j][0] = u[0][cnt]
                    self.r[i][j][1] = u[1][cnt]
                    cnt += 1

            iteration_cnt += 1
            residual = EllipticGrid2D.calc_diff(self.r, old_grid)
            print("{}: {}".format(iteration_cnt, residual))


class EllipticGrid3D(object):
    def __init__(self, grid):
        self.r = np.copy(grid)
        self.d_xi = 1.0
        self.d_eta = 1.0
        self.d_zeta = 1.0

    def r_xi(self, i, j, k):
        return 0.5 * (self.r[i + 1][j][k] - self.r[i - 1][j][k]) / self.d_xi

    def r_eta(self, i, j, k):
        return 0.5 * (self.r[i][j + 1][k] - self.r[i][j - 1][k]) / self.d_eta

    def r_zeta(self, i, j, k):
        return 0.5 * (self.r[i][j][k + 1] - self.r[i][j][k - 1]) / self.d_zeta

    def r_xi2(self, i, j, k):
        return (self.r[i - 1][j][k] - 2 * self.r[i][j][k] + self.r[i + 1][j][k]) / math.pow(self.d_xi, 2)

    def r_eta2(self, i, j, k):
        return (self.r[i][j - 1][k] - 2 * self.r[i][j][k] + self.r[i][j + 1][k]) / math.pow(self.d_eta, 2)

    def r_zeta2(self, i, j, k):
        return (self.r[i][j][k - 1] - 2 * self.r[i][j][k] + self.r[i][j][k + 1]) / math.pow(self.d_zeta, 2)

    def r_xi_eta(self, i, j, k):
        return 0.25 * (self.r[i + 1][j + 1][k] - self.r[i + 1][j - 1][k] - self.r[i - 1][j + 1][k] + self.r[i - 1][j - 1][k]) / (self.d_xi * self.d_eta)

    def r_eta_zeta(self, i, j, k):
        return 0.25 * (self.r[i][j + 1][k + 1] - self.r[i][j - 1][k + 1] - self.r[i][j + 1][k - 1] + self.r[i][j - 1][k - 1]) / (self.d_eta * self.d_zeta)

    def r_zeta_xi(self, i, j, k):
        return 0.25 * (self.r[i + 1][j][k + 1] - self.r[i - 1][j][k + 1] - self.r[i + 1][j][k - 1] + self.r[i - 1][j][k - 1]) / (self.d_zeta * self.d_xi)

    def alpha1(self, i, j, k):
        pass

    def alpha2(self, i, j, k):
        pass

    def alpha3(self, i, j, k):
        pass

    def beta12(self, i, j, k):
        pass

    def beta23(self, i, j, k):
        pass

    def beta31(self, i, j, k):
        pass

    def jacobi(self, i, j, k):
        pass


class EllipticGridTester(unittest.TestCase):
    def test_2d_laplace(self):
        pass

    def test_2d_tm(self):
        pass

    def test_3d_laplace(self):
        pass

    def test_3d_tm(self):
        pass


if __name__ == '__main__':
    unittest.main()
