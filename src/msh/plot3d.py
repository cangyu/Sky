import numpy as np


class PLOT3D(object):
    def __init__(self):
        """
        多块结构网格
        """

        self.blk_list = []

    @property
    def blk_num(self):
        """
        所含块数量
        """

        return len(self.blk_list)

    def add_block(self, blk):
        """
        增加Block
        :param blk: 单块网格
        :type blk: PLOT3D_Block
        :return: None
        """

        self.blk_list.append(blk)

    def write(self, fn, with_iblank=False):
        """
        输出多块网格到文件
        :param fn: 输出文件名
        :param with_iblank: 是否包含IBLANK信息
        :return: None
        """

        fout = open(fn, 'w')
        fout.write("{}".format(self.blk_num))
        for blk in self.blk_list:
            fout.write("\n{} {} {}".format(blk.I, blk.J, blk.K))
        for blk in self.blk_list:
            blk.write(fout, with_iblank)
        fout.close()


"""
IBLANK Values:
 0 : 计算域之外
 1 : 正常点
 2 : 固面边界
-n : 与第n块网格相邻
"""


class PLOT3D_Block(object):
    def __init__(self, pts, idx=0):
        """
        单块结构网格
        :param pts: 所有网格点坐标，下标从左到右依次循环I, J, K, 每个元素包含(X, Y, Z, IBLANK)
        :param idx: 块序号, ranging from 1
        """

        self.index = idx
        self.data = np.copy(pts)

    @property
    def I(self):
        """
        :return X方向节点数量
        :rtype int
        """

        return self.data.shape[0]

    @property
    def J(self):
        """
        :return Y方向节点数量
        :rtype int
        """

        return self.data.shape[1]

    @property
    def K(self):
        """
        :return Z方向节点数量
        :rtype int
        """

        return self.data.shape[2]

    @property
    def pnt_num(self):
        """
        :return 该Block中point的数量
        :rtype int
        """

        return self.I * self.J * self.K

    @property
    def cell_num(self):
        """
        :return 该Block中cell的数量
        :rtype int
        """

        return (self.I - 1) * (self.J - 1) if self.K == 1 else (self.I - 1) * (self.J - 1) * (self.K - 1)

    @property
    def face_num(self):
        """
        :return 该Block中face的数量
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
                    ret[t] = np.copy(self.data[i][j][k][:3])
                    t += 1

        return ret

    def write(self, fout, with_iblank):
        """
        输出网格到流
        :param fout:输出流
        :param with_iblank: 是否包含IBLANK信息
        :return: None
        """

        for d in range(3):
            for k in range(self.K):
                for j in range(self.J):
                    for i in range(self.I):
                        fout.write("{}{}".format('\n' if i == 0 else ' ', self.data[i][j][k][d]))

        if with_iblank:
            for k in range(self.K):
                for j in range(self.J):
                    for i in range(self.I):
                        fout.write("{}{}".format('\n' if i == 0 else ' ', int(self.data[i][j][k][-1])))

    def set_iblank(self, i, j, k, t):
        """
        设置网格点的IBLANK信息
        :param i: X方向下标
        :param j: Y方向下标
        :param k: Z方向下标
        :param t: IBLANK Value
        :type t: int
        :return: None
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
    def build_from_2d(cls, pts, close=True):
        """
        从二维网格点数组构建PLOT3D块
        :param pts: 二维网格点数组, 下标从左到右依次循环I, J, 每个元素包含(X, Y, Z) / (X, Y)
        :param close: 是否封闭边界
        :return: PLOT3D_Block
        """

        I, J, Dim = pts.shape
        p3d = np.zeros((I, J, 1, 4))
        for k in range(1):
            for j in range(J):
                for i in range(I):
                    for d in range(Dim):
                        p3d[i][j][k][d] = pts[i][j][d]
                    p3d[i][j][k][-1] = 1  # 默认均设为正常点

        ret = cls(p3d)
        if close:
            ret.set_boundary_iblank(2)
        return ret

    @classmethod
    def build_from_3d(cls, pts, close=True):
        """
        从三维网格点数组构建PLOT3D块
        :param pts: 三维网格点数组, 下标从左到右依次循环I, J, K, 每个元素包含(X, Y, Z)
        :param close: 是否封闭边界
        :return: PLOT3D_Block
        """

        I, J, K, Dim = pts.shape
        p3d = np.zeros((I, J, K, 4))
        for k in range(K):
            for j in range(J):
                for i in range(I):
                    for d in range(Dim):
                        p3d[i][j][k][d] = pts[i][j][k][d]
                    p3d[i][j][k][-1] = 1  # 默认均设为正常点

        ret = cls(p3d)
        if close:
            ret.set_boundary_iblank(2)
        return ret
