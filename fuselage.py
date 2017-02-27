from abc import ABCMeta, abstractmethod
import matplotlib.pyplot as plt
import numpy as np
import math, os, sys
from stl import mesh
from PyQt5.QtWidgets import *


# 机身抽象定义
class Fuselage():
    __metaclass__ = ABCMeta

    def __init__(self, _filename):
        # real description, 右手系, 机头鼻尖为(0, 0, 0)
        self.filename = _filename  # 当前机身的文件名
        self.x = []  # 各个剖面沿x方向的坐标
        self.y_up = []  # 机身在xy平面内投影，上方曲线的y坐标
        self.y_down = []  # 机身在xy平面内投影，下方曲线的y坐标
        self.y_mid = []  # 机身在xy平面内投影，中线的y坐标

        # intrinsic description
        self.L = 57  # 全长
        self.D = 5.97  # 中部直径
        self.Theta_fc = 14  # 擦地角

        self.r0 = 0.15  # 头部长度比例
        self.r2 = 0.26  # 尾部长度比例
        self.r1 = 1.0 - self.r0 - self.r2  # 中部长度比例

        # derived description
        self.L0 = self.L * self.r0  # 头部长度
        self.L1 = self.L * self.r1  # 中部长度
        self.L2 = self.L * self.r2  # 尾部长度

        self.lambda_total = self.L / self.D  # 全机长径比
        self.lambda_front = self.L0 / self.D  # 头部长径比
        self.lambda_mid = self.L1 / self.D  # 中部长径比
        self.lambda_tail = self.L2 / self.D  # 尾部长径比

        self.HeadingDirCapacity = 0.125  # 航向机身容量参数
        self.PitchDirCapacity = 1.25  # 纵向机身容量参数

    @abstractmethod
    def calc_y_up(self, x):
        pass

    @abstractmethod
    def calc_y_down(self, x):
        pass

    def calc_profile(self, _step):
        # 生成x方向离散点
        x_front = np.arange(0, self.L0, _step)
        x_mid = np.arange(self.L0, self.L0 + self.L1, _step)
        x_tail = np.arange(self.L0 + self.L1, self.L, _step)

        self.x.clear()

        for i in range(0, len(x_front)):
            self.x.append(x_front[i])

        for i in range(0, len(x_mid)):
            self.x.append(x_mid[i])

        for i in range(0, len(x_tail)):
            self.x.append(x_tail[i])

        self.x.append(self.L)

        # 计算每点对应的上下两个y坐标
        self.y_up.clear()
        self.y_down.clear()
        for i in range(0, len(self.x)):
            self.y_up.append(self.calc_y_up(self.x[i]))
            self.y_down.append(self.calc_y_down(self.x[i]))

    def calc_mid_line(self):
        self.y_mid.clear()
        assert len(self.x) > 0

        assert len(self.x) == len(self.y_up)
        assert len(self.x) == len(self.y_down)

        for i in range(0, len(self.x)):
            self.y_mid.append(float(self.y_up[i] + self.y_down[i]) / 2)

    def show_profile(self):
        self.calc_mid_line()
        plt.plot(self.x, self.y_up, label="UP")
        plt.plot(self.x, self.y_down, label="DOWN")
        plt.plot(self.x, self.y_mid, label="MID")
        plt.gca().set_aspect(1)
        plt.show()

    def generate(self, _delta_len=0.05, _delta_theta=10):
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
        self.show_profile()


# 剖面头尾为1/4椭圆组成，不带擦地角
class SimpleFuselage(Fuselage):
    def __init__(self, _filename):
        Fuselage.__init__(self, _filename)

        # intrinsic description
        self.hr0 = 0.8  # 机头上部高度与直径之比
        self.hr2 = 0.18  # 机尾上部高度与直径之比

        # derived description
        self.h0 = self.D * self.hr0  # 机头上部高度
        self.h1 = self.D - self.h0  # 机头下部高度
        self.h2 = self.D * self.hr2  # 机尾上部高度
        self.h3 = self.D - self.h2  # 机尾下部高度
        self.dh = self.h0 - self.h2  # 尾缘顶点与前缘顶点之间的高度差

        # widget
        self.hr0_label = QLabel('机头上层占直径比例:')
        self.hr0_dsb = QDoubleSpinBox()
        self.hr0_dsb.setRange(0.0, 1.0)
        self.hr0_dsb.setSingleStep(0.01)
        self.hr0_dsb.setValue(self.hr0)

        self.hr2_label = QLabel('机尾上层占直径比例:')
        self.hr2_dsb = QDoubleSpinBox()
        self.hr2_dsb = QDoubleSpinBox()
        self.hr2_dsb.setRange(0.0, 1.0)
        self.hr2_dsb.setSingleStep(0.01)
        self.hr2_dsb.setValue(self.hr2)

    def init_widget(self):
        Fuselage.init_widget(self)

        self.intrinsic_param_layout.addWidget(self.hr0_label, self.intrinsic_param_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.hr0_dsb, self.intrinsic_param_cnt, 1)
        self.intrinsic_param_cnt += 1

        self.intrinsic_param_layout.addWidget(self.hr2_label, self.intrinsic_param_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.hr2_dsb, self.intrinsic_param_cnt, 1)
        self.intrinsic_param_cnt += 1

    def update_derived_param(self):
        Fuselage.update_derived_param(self)

        self.hr0 = self.hr0_dsb.value()
        self.hr2 = self.hr2_dsb.value()

        self.h0 = self.D * self.hr0
        self.h1 = self.D - self.h0
        self.h2 = self.D * self.hr2
        self.h3 = self.D - self.h2
        self.dh = self.h0 - self.h2

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


# 剖面头为两个1/4椭圆组成，带擦地角
class TraditionalFuselage(Fuselage):
    def __init__(self, _filename):
        Fuselage.__init__(self, _filename)

        # intrinsic description
        self.Theta_fc = 14.0
        self.hr0 = 0.84  # 机头上部高度与直径之比
        self.de = 0.7  # 机尾圆圈直径

        # derived description
        self.h0 = self.D * self.hr0  # 机头上部高度
        self.h1 = self.D - self.h0  # 机头下部高度

        self.dh = math.tan(math.radians(self.Theta_fc)) * self.L2 + self.de - self.h1
        self.h2 = self.D - self.h1 - self.dh

        # widget
        self.hr0_label = QLabel('机头上层占直径比例:')
        self.hr0_dsb = QDoubleSpinBox()
        self.hr0_dsb.setRange(0.0, 1.0)
        self.hr0_dsb.setSingleStep(0.01)
        self.hr0_dsb.setValue(self.hr0)

        self.de_label = QLabel('机尾圆圈直径(m):')
        self.de_dsb = QDoubleSpinBox()
        self.de_dsb = QDoubleSpinBox()
        self.de_dsb.setRange(0.0, 1.0)
        self.de_dsb.setSingleStep(0.01)
        self.de_dsb.setValue(self.de)

    def init_widget(self):
        Fuselage.init_widget(self)

        self.intrinsic_param_layout.addWidget(self.hr0_label, self.intrinsic_param_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.hr0_dsb, self.intrinsic_param_cnt, 1)
        self.intrinsic_param_cnt += 1

        self.intrinsic_param_layout.addWidget(self.de_label, self.intrinsic_param_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.de_dsb, self.intrinsic_param_cnt, 1)
        self.intrinsic_param_cnt += 1

    def update_derived_param(self):
        Fuselage.update_derived_param(self)

        self.hr0 = self.hr0_dsb.value()
        self.de = self.de_dsb.value()

        self.h0 = self.D * self.hr0
        self.h1 = self.D - self.h0

        self.dh = math.tan(math.radians(self.Theta_fc)) * self.L2 + self.de - self.h1
        self.h2 = self.D - self.h1 - self.dh

    def calc_y_up(self, x):
        if x <= 0:
            return 0
        elif x <= self.L0:
            return self.h0 * math.sqrt(1 - math.pow((x - self.L0) / self.L0, 2))
        elif x <= self.L0 + self.L1:
            return self.h0
        elif x <= self.L:
            return self.dh + self.h2 * math.sqrt(1 - math.pow((x - (self.L0 + self.L1)) / self.L2, 2))
        else:
            return self.dh

    def calc_y_down(self, x):
        if x <= 0:
            return 0
        elif x <= self.L0:
            return -self.h1 * math.sqrt(1 - math.pow((x - self.L0) / self.L0, 2))
        elif x <= self.L0 + self.L1:
            return -self.h1
        elif x <= self.L:
            return -self.h1 + math.tan(math.radians(self.Theta_fc)) * (x - (self.L0 + self.L1))
        else:
            return 0
