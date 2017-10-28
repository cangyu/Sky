import numpy as np
import math
import unittest

"""
Implementation of the Plot3D standard.
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
    def I(self):
        """
        :return Dimension in the 1st direction.
        :rtype int
        """

        return self.data.shape[0]

    @property
    def J(self):
        """
        :return Dimension in the 2nd direction.
        :rtype int
        """

        return self.data.shape[1]

    @property
    def K(self):
        """
        :return Dimension in the 3rd direction.
        :rtype int
        """

        return self.data.shape[2]

    def __repr__(self):
        return "{} x {}".format(self.I, self.J) if self.K == 1 else "{} x {} x {}".format(self.I, self.J, self.K)

    @property
    def pnt_num(self):
        """
        :return Number of points.
        :rtype int
        """

        return self.I * self.J * self.K

    @property
    def cell_num(self):
        """
        :return Number of cells.
        :rtype int
        """

        return (self.I - 1) * (self.J - 1) if self.K == 1 else (self.I - 1) * (self.J - 1) * (self.K - 1)

    @property
    def face_num(self):
        """
        :return Number of faces.
        :rtype int
        """

        return 3 * self.I * self.J * self.K - (self.I * self.J + self.J * self.K + self.K * self.I)

    @property
    def all_pnt(self):
        t = 0
        ret = np.empty((self.pnt_num, 3), float)
        for k in range(self.K):
            for j in range(self.J):
                for i in range(self.I):
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

        for d in range(3):
            for k in range(self.K):
                for j in range(self.J):
                    for i in range(self.I):
                        f_out.write("{}{}".format('\n' if i == 0 else ' ', self.data[i][j][k][d]))

        if with_iblank:
            for k in range(self.K):
                for j in range(self.J):
                    for i in range(self.I):
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

        if self.K == 1:
            for i in range(self.I):
                for j in range(self.J):
                    if i in (0, self.I - 1) or j in (0, self.J - 1):
                        self.set_iblank(i, j, 0, t)
        else:
            for i in range(self.I):
                for j in range(self.J):
                    for k in range(self.K):
                        if i in (0, self.I - 1) or j in (0, self.J - 1) or k in (0, self.K - 1):
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


if __name__ == '__main__':
    unittest.main()
