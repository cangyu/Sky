import math
import numpy as np
import matplotlib.pyplot as plt
from stl import mesh


# 简单机身
class Fuselage(object):
    def __init__(self, _filename):
        # real description
        self.filename = _filename
        self.x = []
        self.y_up = []
        self.y_down = []
        self.y_mid = []

        # intrinsic description
        self.L = 6380  # 全长
        self.D = 700  # 中部直径
        self.Theta_fc = 30  # 擦地角
        self.r0 = 0.18 #头部长度比例
        self.r2 = 0.23 #尾部长度比例
        self.r1 = 1.0 - self.r0 - self.r2 #中部长度比例

        self.L0 = self.L * self.r0 #头部长度
        self.L1 = self.L * self.r1 #中部长度
        self.L2 = self.L * self.r2 #尾部长度


    def set_parameter(self, _l, _d, _tfc, _r0, _r2, _hr0, _hr2):
        '''
        外形参数,右手系,机头鼻尖为(0,0,0)
        :param _l: 全机身长度
        :param _d: 中部等截面部分直径
        :param _tfc: 擦地角
        :param _r0: 机头长度与全机长度之比
        :param _r2: 机尾长度与全机长度之比
        :param _hr0: 机头上部高度与直径之比
        :param _hr2: 机尾上部高度与直径之比
        :return:
        '''

        assert _l > 0
        self.L = _l

        assert _d > 0
        self.D = _d

        assert _tfc > 0
        self.Theta_fc = _tfc

        assert 0 < _r0 and _r0 < 1.0
        self.r0 = _r0

        assert 0 < _r2 and _r2 < 1.0
        self.r2 = _r2

        assert _r0 + _r2 < 1.0
        self.r1 = 1.0 - _r0 - _r2

        self.L0 = self.L * self.r0
        self.L1 = self.L * self.r1
        self.L2 = self.L * self.r2

        assert 0 < _hr0 and _hr0 < 1.0
        self.hr0 = _hr0

        assert 0 < _hr2 and _hr2 < 1.0
        self.hr2 = _hr2

        self.h0 = self.D * _hr0
        self.h1 = self.D - self.h0

        self.h2 = self.D * _hr2
        self.h3 = self.D - self.h2

        self.dh = self.h0 - self.h2

    def calc_relevant_parameter(self):
        self.lambda_total = self.L / self.D
        self.lambda_front = self.L0 / self.D
        self.lambda_mid = self.L1 / self.D
        self.lambda_tail = self.L2 / self.D

        return self.lambda_total, self.lambda_front, self.lambda_mid, self.lambda_tail

    def calc_y_up(self, x):
        if x < 0:
            return 0
        elif x < self.L0:
            return self.h0 * math.sqrt(1 - math.pow((x - self.L0) / self.L0, 2))
        elif x < self.L0 + self.L1:
            return self.h0
        elif x < self.L:
            return self.dh + self.h2 * math.sqrt(1 - math.pow((x - (self.L0 + self.L1)) / self.L2, 2))
        else:
            return self.dh

    def calc_y_down(self, x):
        if x < 0:
            return 0
        elif x < self.L0:
            return -self.h1 * math.sqrt(1 - math.pow((x - self.L0) / self.L0, 2))
        elif x < self.L0 + self.L1:
            return -self.h1
        elif x < self.L:
            return self.dh - self.h3 * math.sqrt(1 - math.pow((x - (self.L0 + self.L1)) / self.L2, 2))
        else:
            return self.dh

    def calc_profile(self, _step):
        # x方向离散点
        x_front = np.arange(0, self.L0, _step)
        x_mid = np.arange(self.L0, self.L0 + self.L1, _step)
        x_tail = np.arange(self.L0 + self.L1, self.L, _step)

        for i in range(0, len(x_front)):
            self.x.append(x_front[i])

        for i in range(0, len(x_mid)):
            self.x.append(x_mid[i])

        for i in range(0, len(x_tail)):
            self.x.append(x_tail[i])

        self.x.append(self.L)

        for i in range(0, len(self.x)):
            self.y_up.append(self.calc_y_up(self.x[i]))
            self.y_down.append(self.calc_y_down(self.x[i]))

    def calc_mid_line(self):
        self.y_mid.clear()
        for i in range(0, len(self.x)):
            self.y_mid.append((self.y_up[i] + self.y_down[i]) / 2)

    def show_profile(self):
        self.calc_mid_line()
        plt.plot(self.x, self.y_up, label="UP")
        plt.plot(self.x, self.y_down, label="DOWN")
        plt.plot(self.x, self.y_mid, label="MID")
        plt.gca().set_aspect(1)
        plt.show()

    def generate(self, _delta_len=5, _delta_theta=10):
        fuselage_surf = []

        self.calc_profile(_delta_len)

        for i in range(0, len(self.x) - 1):
            rl = (self.y_up[i] - self.y_down[i]) / 2
            rr = (self.y_up[i + 1] - self.y_down[i + 1]) / 2

            theta = 0
            delta_theta = _delta_theta
            while theta < 360:
                pl0 = [self.x[i],
                       rl * math.cos(math.radians(theta)),
                       rl * math.sin(math.radians(theta))]
                pl1 = [self.x[i],
                       rl * math.cos(math.radians(theta + delta_theta)),
                       rl * math.sin(math.radians(theta + delta_theta))]
                pr0 = [self.x[i + 1],
                       rr * math.cos(math.radians(theta)),
                       rr * math.sin(math.radians(theta))]
                pr1 = [self.x[i + 1],
                       rr * math.cos(math.radians(theta + delta_theta)),
                       rr * math.sin(math.radians(theta + delta_theta))]

                fuselage_surf.append([pl0, pl1, pr0])
                fuselage_surf.append([pr0, pr1, pl1])

                theta = theta + delta_theta

        fuselage = mesh.Mesh(np.zeros(len(fuselage_surf), dtype=mesh.Mesh.dtype))
        for i in range(0, len(fuselage_surf)):
            fuselage.vectors[i] = fuselage_surf[i]

        fuselage.save(self.filename)

        return fuselage


# 剖面头尾为1/4椭圆组成
class SimpleFuselage(Fuselage):
    def __init__(self, _filename):
        Fuselage.__init__(self, _filename)

        self.hr0 = 0.8
        self.hr2 = 0.18

        self.h0 = self.D * self.hr0
        self.h1 = self.D - self.h0

        self.h2 = self.D * self.hr2
        self.h3 = self.D - self.h2

        self.dh = self.h0 - self.h2

