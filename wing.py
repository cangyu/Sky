import matplotlib.pyplot as plt
import numpy as np
import math, os
from stl import mesh


# 2维翼型描述, 弦长为单位1
class Airfoil(object):
    def __init__(self, _file):
        self.x = []
        self.y_up = []
        self.y_down = []

        airfoil = open(_file)
        for line in airfoil:
            (_x, _y_up, _y_down) = line.split()
            self.x.append(float(_x))
            self.y_up.append(float(_y_up))
            self.y_down.append(float(_y_down))
        airfoil.close()

    def show(self):
        plt.plot(self.x, self.y_up, label="UP")
        plt.plot(self.x, self.y_down, label="DOWN")
        plt.show()


# 组成3维机翼的真实剖面
class Section(object):
    def __init__(self, _airfoil, _z, _x_front, _x_tail, _dy, _theta):
        self.airfoil = _airfoil  # 剖面的翼型
        self.z = _z  # 剖面与根部的距离
        self.x_front = _x_front  # 前缘点在x方向上的位置
        self.x_tail = _x_tail  # 后缘点在x方向上的位置
        self.dy = _dy  # 由于机翼上反角导致的在y方向上引起的偏移
        self.theta = _theta  # 绕前缘的扭转角

        self.n = len(_airfoil.x)
        self.chord = np.zeros((self.n, 3))  # 剖面弦线上的离散点
        self.up_pts = np.zeros((self.n, 3))  # 剖面上表面的点
        self.down_pts = np.zeros((self.n, 3))  # 剖面下表面的点

    # 计算描述剖面形状的空间坐标点
    def calc_discrete_pts(self):

        # 弦长
        chord_len = self.x_tail - self.x_front

        # 计算正常放置时上、下表面上每点的坐标
        for i in range(0, self.n):
            # x方向
            self.up_pts[i][0] = self.x_front + chord_len * self.airfoil.x[i]
            self.down_pts[i][0] = self.up_pts[i][0]
            # y方向
            self.up_pts[i][1] = chord_len * self.airfoil.y_up[i]
            self.down_pts[i][1] = chord_len * self.airfoil.y_down[i]
            # z方向
            self.up_pts[i][2] = self.z
            self.down_pts[i][2] = self.z

        # 计算由于绕前缘扭转带来的改变
        rotate = complex(math.cos(math.radians(self.theta)), math.sin(math.radians(self.theta)))
        for i in range(1, self.n):
            t = complex(self.up_pts[i][0] - self.up_pts[0][0], self.up_pts[i][1] - self.up_pts[0][1]) * rotate
            self.up_pts[i][0] = self.up_pts[0][0] + t.real
            self.up_pts[i][1] = self.up_pts[0][1] + t.imag
        for i in range(1, self.n):
            t = complex(self.down_pts[i][0] - self.down_pts[0][0], self.down_pts[i][1] - self.down_pts[0][1]) * rotate
            self.down_pts[i][0] = self.down_pts[0][0] + t.real
            self.down_pts[i][1] = self.down_pts[0][1] + t.imag

        # 加上由于上反角导致的y方向上的增量
        for i in range(0, self.n):
            self.up_pts[i][1] += self.dy
            self.down_pts[i][1] += self.dy

        # 计算弦线上每点坐标，3维
        for i in range(0, self.n):
            self.chord[i][0] = 0.5 * (self.up_pts[i][0] + self.down_pts[i][0])
            self.chord[i][1] = 0.5 * (self.up_pts[i][1] + self.down_pts[i][1])
            self.chord[i][2] = 0.5 * (self.up_pts[i][2] + self.down_pts[i][2])

    def generate_face(self):
        face = []
        face.append([self.up_pts[0], self.up_pts[1], self.chord[1]])
        face.append([self.down_pts[0], self.down_pts[1], self.chord[1]])
        for i in range(1, self.n - 1):
            face.append([self.chord[i], self.up_pts[i], self.up_pts[i + 1]])
            face.append([self.chord[i], self.down_pts[i], self.down_pts[i + 1]])
            face.append([self.chord[i], self.chord[i + 1], self.up_pts[i + 1]])
            face.append([self.chord[i], self.chord[i + 1], self.down_pts[i + 1]])

        return face

    def generate_surf(cls, start, end):
        '''
        生成两个剖面之间的曲面，两个剖面要有相同的样本点数量n
        :param start: 起始剖面
        :param end: 结束剖面
        :return: 一个数组，包含了构成曲面的所有三角形
        '''
        assert start.n == end.n

        surf = []
        n = start.n
        # 上表面
        for i in range(0, n - 1):
            surf.append([start.up_pts[i], start.up_pts[i + 1], end.up_pts[i]])
            surf.append([end.up_pts[i], end.up_pts[i + 1], start.up_pts[i + 1]])

        # 下表面
        for i in range(0, n - 1):
            surf.append([start.down_pts[i], start.down_pts[i + 1], end.down_pts[i]])
            surf.append([end.down_pts[i], end.down_pts[i + 1], start.down_pts[i + 1]])

        # 后缘垂面
        surf.append([start.up_pts[n - 1], start.down_pts[n - 1], end.up_pts[n - 1]])
        surf.append([end.up_pts[n - 1], end.down_pts[n - 1], start.down_pts[n - 1]])

        return surf


# 简单机翼, 由一个个剖面组成
class Wing(object):

    def __init__(self, _filename):

        # real description
        self.filename = _filename  # 文件名
        self.SectionNum = 12  # 剖面数量
        self.section_list = []  # 剖面列表
        self.airfoil_list = []  # 翼型列表
        self.z_list = []  # 各个剖面在z轴的位置列表
        self.x_front_list = []  # 各个剖面的前缘点在x轴上的位置
        self.x_tail_list = []  # 各个剖面的尾缘点在x轴上的位置
        self.dy_list = []  # 各个剖面在y轴上偏移量
        self.theta_list = []  # 各个剖面绕前缘点扭转角度

        # intrinsic description
        self.Airfoil = ''  # 翼型
        self.Span = 30  # 展长
        self.C_root = 8.0  # 翼根弦长
        self.C_tip = 0.8  # 翼尖弦长
        self.SweepBack = 32.2  # 后掠角
        self.Dihedral = 5  # 上反角
        self.Twist = -4  # 绕前缘扭转角

        # installation description
        self.X_25 = 17.22  # 1/4弦长位置
        self.dY = -1.0  # Y方向偏移量
        self.dZ = 0  # Z方向偏移量

        # derived description
        self.S = 0  # 参考机翼面积
        self.AR = 0  # 展弦比
        self.A_25 = 0  # 1/4弦线后掠角
        self.TaperRatio = 0  # 梯形比
        self.MAC = 0  # 平均气动弦长

    def set_param(self, _ui):
        self.Airfoil = _ui.Airfoil
        self.Span = _ui.Span
        self.C_root = _ui.C_root
        self.C_tip = _ui.C_tip
        self.SweepBack = _ui.SweepBack
        self.Dihedral = _ui.Dihedral
        self.Twist = _ui.Twist
        self.X_25 = _ui.X_25
        self.dY = _ui.dY
        self.dZ = _ui.dZ

    def set_real_description(self, _airfoil, _z, _x_front, _x_tail, _dy, _theta):
        self.SectionNum = len(_z)
        self.airfoil_list = _airfoil
        self.z_list = _z
        self.x_front_list = _x_front
        self.x_tail_list = _x_tail
        self.dy_list = _dy
        self.theta_list = _theta

    def calc_sections(self):
        self.section_list.clear()

        for i in range(0, self.SectionNum):
            self.section_list.append(
                Section(self.airfoil_list[i], self.z_list[i], self.x_front_list[i], self.x_tail_list[i],
                        self.dy_list[i], self.theta_list[i]))

        for i in range(0, len(self.z_list)):
            self.section_list[i].calc_discrete_pts()

    # 目前只生成全部参数线性分布的简单机翼
    def calc_real_desc(self):

        # 剖面翼型
        _airfoil = Airfoil(self.Airfoil)
        airfoil_dist = []
        for i in range(0, self.SectionNum):
            airfoil_dist.append(_airfoil)

        # 剖面分布
        z_dist = np.linspace(0, self.Span, self.SectionNum)

        # 前缘点
        x_front_dist = np.linspace(0, self.Span * math.tan(math.radians(self.SweepBack)), self.SectionNum)

        # 后缘点
        x_tail_dist = []
        length = np.linspace(self.C_root, self.C_tip, self.SectionNum)
        for i in range(0, self.SectionNum):
            x_tail_dist.append(x_front_dist[i] + length[i])

        # 上反
        dy_dist = np.linspace(0, self.Span * math.tan(math.radians(self.Dihedral)), self.SectionNum)

        # 扭转
        theta_dist = []
        for i in range(0, self.SectionNum):
            theta_dist.append(self.Twist)

        # Wing generation
        self.set_real_description(airfoil_dist, z_dist, x_front_dist, x_tail_dist, dy_dist, theta_dist)
        self.calc_sections()

    def show_profile(self):
        x_25=[]

        for i in range(0, len(self.z_list)):
            x_25.append(float(0.75*self.x_front_list[i]+0.25*self.x_tail_list[i]))

        plt.plot(self.z_list, self.x_front_list, label="front")
        plt.plot(self.z_list, x_25, label="mid")
        plt.plot(self.z_list, self.x_tail_list, label="tail")
        plt.gca().set_aspect(1)
        plt.show()

    def generate(self):
        self.calc_real_desc()
        wing_surf = []

        # 上下表面与尾缘
        for i in range(0, len(self.section_list) - 1):
            for face in self.section_list[i].generate_surf(self.section_list[i], self.section_list[i + 1]):
                wing_surf.append(face)

        # 根部与翼梢
        for face in self.section_list[0].generate_face():
            wing_surf.append(face)
        for face in self.section_list[-1].generate_face():
            wing_surf.append(face)

        # 生成STL格式文件
        wing = mesh.Mesh(np.zeros(len(wing_surf), dtype=mesh.Mesh.dtype))
        for i in range(0, len(wing_surf)):
            wing.vectors[i] = wing_surf[i]

        wing.save(self.filename)
        self.show_profile()


# 垂尾
class VerticalStabilizer(Wing):
    def __init__(self, _filename):
        Wing.__init__(self, _filename)

        # real description
        self.SectionNum = 5

        # intrinsic description
        self.Span = 9.0
        self.C_root = 6.6
        self.C_tip = 2.6
        self.SweepBack = 40
        self.Dihedral = 0
        self.Twist = 0

        # installation param
        self.X_25 = 50.0
        self.dY = 0
        self.dZ = 0

        # derived description
        self.V_v = 0.079  # 垂尾容量


# 平尾
class HorizontalStabilizer(Wing):
    def __init__(self, _filename):
        Wing.__init__(self, _filename)

        # real description
        self.SectionNum = 4

        # intrinsic description
        self.Span = 10.2
        self.C_root = 4.5
        self.C_tip = 0.9
        self.SweepBack = 25
        self.Dihedral = 0
        self.Twist = 0

        # installation description
        self.X_25 = 50.0
        self.dY = 0
        self.dZ = 0

        # derived description
        self.V_h = 0
