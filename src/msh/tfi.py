import numpy as np


class Linear_TFI_2D(object):
    def __init__(self, c1, c2, c3, c4):
        """
        2维无限插值(Linear)
        :param c1: 平行于x轴方向的第1条曲线，调用得到3维坐标点
        :param c2: 平行于y轴方向的第1条曲线，调用得到3维坐标点
        :param c3: 平行于x轴方向的第2条曲线，调用得到3维坐标点
        :param c4: 平行于y轴方向的第2条曲线，调用得到3维坐标点
        """

        self.grid = None

        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4

        self.P12 = c1(0)  # c1与c2交点
        self.P14 = c1(1)  # c1与c4交点
        self.P23 = c3(0)  # c2与c3交点
        self.P34 = c3(1)  # c3与c4交点

        self.U = lambda u, v: (1 - u) * self.c2(v) + u * self.c4(v)
        self.V = lambda u, v: (1 - v) * self.c1(u) + v * self.c3(u)
        self.UV = lambda u, v: (1 - u) * (1 - v) * self.P12 + u * v * self.P34 + (1 - u) * v * self.P23 + (1 - v) * u * self.P14

    def __call__(self, u, v):
        """
        曲面在(u,v)处的坐标
        """

        return self.U(u, v) + self.V(u, v) - self.UV(u, v)

    def calc_grid(self, pu, pv):
        """
        根据网格点的参数分布计算对应的坐标
        :param pu: 所有网格点的U方向参数值，2维，IxJ 个网格点
        :param pv: 所有网格点的V方向参数值，2维，IxJ 个网格点
        :return: 所有网格点的坐标
        """

        if len(pu.shape) != 2:
            raise AssertionError("Invalid dimension!")
        if pu.shape != pv.shape:
            raise AssertionError("U V shape don't match!")

        d1, d2 = pu.shape
        self.grid = np.empty((d1, d2, 3), float)

        for i in range(d1):
            for j in range(d2):
                self.grid[i][j] = self.__call__(pu[i][j], pv[i][j])

        return self.grid


class Linear_TFI_3D(object):
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

        self.grid = None

        self.s1 = s1
        self.s2 = s2
        self.s3 = s3
        self.s4 = s4
        self.s5 = s5
        self.s6 = s6

        self.c35 = lambda u: s5(u, 0)  # intersection of s5, s3
        self.c15 = lambda v: s5(0, v)  # intersection of s5, s1
        self.c45 = lambda u: s5(u, 1)  # intersection of s5, s4
        self.c25 = lambda v: s5(1, v)  # intersection of s5, s2
        self.c36 = lambda u: s6(u, 0)  # intersection of s6, s3
        self.c16 = lambda v: s6(0, v)  # intersection of s6, s1
        self.c46 = lambda u: s6(u, 1)  # intersection of s6, s4
        self.c26 = lambda v: s6(1, v)  # intersection of s6, s2
        self.c13 = lambda w: s3(w, 0)  # intersection of s1, s3
        self.c23 = lambda w: s2(0, w)  # intersection of s3, s2
        self.c24 = lambda w: s4(w, 1)  # intersection of s2, s4
        self.c14 = lambda w: s1(1, w)  # intersection of s4, s1

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
        self.UVW = lambda u, v, w: (1 - u) * (1 - v) * (1 - w) * self.p135 + \
                                   (1 - u) * (1 - v) * w * self.p136 + \
                                   (1 - u) * v * (1 - w) * self.p145 + \
                                   (1 - u) * v * w * self.p146 + \
                                   u * (1 - v) * (1 - w) * self.p235 + \
                                   u * (1 - v) * w * self.p236 + \
                                   u * v * (1 - w) * self.p245 + \
                                   u * v * w * self.p246

    def __call__(self, u, v, w):
        """
        在(u,v,w)处的坐标
        """

        return self.U(u, v, w) + self.V(u, v, w) + self.W(u, v, w) - (self.UV(u, v, w) + self.VW(u, v, w) + self.WU(u, v, w)) + self.UVW(u, v, w)

    def calc_grid(self, pu, pv, pw):
        """
        根据网格点的参数分布计算对应的坐标
        :param pu: 所有网格点的U方向参数值，3维，IxJxK 个网格点
        :param pv: 所有网格点的V方向参数值，3维，IxJxK 个网格点
        :param pw: 所有网格点的W方向参数值，3维，IxJxK 个网格点
        :return: 所有网格点的坐标
        """

        if len(pu.shape) != 3:
            raise AssertionError("Invalid dimension!")
        if pu.shape != pv.shape:
            raise AssertionError("U V shape don't match!")
        if pv.shape != pw.shape:
            raise AssertionError("V W shape don't match!")

        d1, d2, d3 = pu.shape
        self.grid = np.empty((d1, d2, d3, 3))

        for i in range(d1):
            for j in range(d2):
                for k in range(d3):
                    self.grid[i][j][k] = self.__call__(pu[i][j][k], pv[i][j][k], pw[i][j][k])

        return self.grid
