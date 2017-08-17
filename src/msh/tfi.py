import numpy as np
from copy import deepcopy


class TFI(object):
    def __init__(self):
        self.grid = None

    def get_grid(self):
        return self.grid


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

        self.c1 = deepcopy(c1)
        self.c2 = deepcopy(c2)
        self.c3 = deepcopy(c3)
        self.c4 = deepcopy(c4)

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
