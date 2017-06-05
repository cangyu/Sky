import numpy as np


class PLOT3D(object):
    def __init__(self):
        """
        多块结构网格
        :param blk: 块列表 
        """

        self.blk_num = 0
        self.blk_list = []

    def add_block(self, blk):
        self.blk_list.append(blk)
        self.blk_num += 1

    def write(self, fn):
        """
        输出多块网格到文件
        :param fn: 输出文件名
        :return: None
        """

        fout = open(fn, 'w')
        fout.write("{}".format(self.blk_num))
        for blk in self.blk_list:
            fout.write("\n{} {} {}".format(blk.I, blk.J, blk.K))
        for blk in self.blk_list:
            blk.write(fout)
        fout.close()


class PLOT3D_Block(object):
    def __init__(self, pts):
        """
        单块结构网格
        :param pts: 所有网格点坐标，下标从左到右依次循环X, Y, Z
        """

        '''X方向节点数量，Y方向节点数量，Z方向节点数量，网格点信息维度'''
        self.I, self.J, self.K, self.Dim = pts.shape
        self.data = np.empty((self.I, self.J, self.K, self.Dim), float, order='F')

        for d in range(self.Dim):
            for k in range(self.K):
                for j in range(self.J):
                    for i in range(self.I):
                        self.data[i][j][k][d] = pts[i][j][k][d]

    def write(self, fout):
        """
        输出网格到流
        :param fout:输出流 
        :return: None
        """

        for d in range(self.Dim):
            for k in range(self.K):
                for j in range(self.J):
                    for i in range(self.I):
                        fout.write("{}{}".format(('\n' if i == 0 else ' '), self.data[i][j][k][d]))

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
        :param t: IBLANK value
        :return: None
        """

        self.data[i][j][k][-1] = t
