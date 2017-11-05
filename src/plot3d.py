import unittest
import math
import numpy as np
from settings import GRID_DBG

"""
Implementation of the Plot3D standard.

Note:
All the units are SI by default.
"""


def report_process(msg):
    if GRID_DBG:
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


if __name__ == '__main__':
    unittest.main()
