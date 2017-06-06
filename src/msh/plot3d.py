import numpy as np


class PLOT3D(object):
    def __init__(self):
        """
        多块结构网格
        """

        self.blk_num = 0
        self.blk_list = []

    def add_block(self, blk):
        self.blk_list.append(blk)
        self.blk_num += 1

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


class PLOT3D_Block(object):
    def __init__(self, pts):
        """
        单块结构网格
        :param pts: 所有网格点坐标，下标从左到右依次循环I, J, K, 
                    每个元素包含(X, Y, Z, IBLANK)
        """

        '''X方向节点数量，Y方向节点数量，Z方向节点数量，网格点信息维度'''
        self.I, self.J, self.K, self.Dim = pts.shape
        self.data = np.empty((self.I, self.J, self.K, self.Dim), float, order='F')

        for d in range(self.Dim):
            for k in range(self.K):
                for j in range(self.J):
                    for i in range(self.I):
                        self.data[i][j][k][d] = pts[i][j][k][d]

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
         0 : 计算域之外
         1 : 正常点
         2 : 固面边界
        -n : 与第n块网格相邻
        :param i: X方向下标
        :param j: Y方向下标
        :param k: Z方向下标
        :param t: IBLANK Value
        :return: None
        """

        self.data[i][j][k][-1] = t

    def set_area_iblank(self, rg0, rg1, rg2, t):
        """
        设置区域内网格点的IBLANK信息
         0 : 计算域之外
         1 : 正常点
         2 : 固面边界
        -n : 与第n块网格相邻
        :param rg0: X方向范围
        :param rg1: Y方向范围
        :param rg2: Z方向范围
        :param t: IBLANK Value
        :return: None
        """

        for i in rg0:
            for j in rg1:
                for k in rg2:
                    self.set_iblank(i, j, k, t)

    def set_boundary_iblank(self, t):
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
        :param pts: 二维网格点数组, 下标从左到右依次循环I, J,
                    每个元素包含(X, Y, Z) / (X, Y)
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

        ret = PLOT3D_Block(p3d)
        if close:
            ret.set_boundary_iblank(2)
        return ret

    @classmethod
    def build_from_3d(cls, pts, close=True):
        """
        从三维网格点数组构建PLOT3D块
        :param pts: 三维网格点数组, 下标从左到右依次循环I, J, K
                    每个元素包含(X, Y, Z)
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

        ret = PLOT3D_Block(p3d)
        if close:
            ret.set_boundary_iblank(2)
        return ret
