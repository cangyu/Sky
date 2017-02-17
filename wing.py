import matplotlib.pyplot as plt
import numpy as np
import math, os
from stl import mesh
from PyQt5.QtWidgets import *


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
    airfoil_list = []  # 翼型列表
    airfoil_dir = './airfoil/'

    @staticmethod
    def update_airfoil_list():
        for f in os.listdir(Wing.airfoil_dir):
            cur_filename = os.path.join(Wing.airfoil_dir, f)
            if os.path.isfile(cur_filename) and cur_filename.endswith('.dat'):
                Wing.airfoil_list.append(f)

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
        self.Span = 34.8  # 展长
        self.C_root = 13.0  # 翼根弦长
        self.C_tip = 2.58  # 翼尖弦长
        self.SweepBack = 30  # 后掠角
        self.Dihedral = 5  # 上反角
        self.Twist = -4  # 绕前缘扭转角

        # installation description
        self.X_25 = 20.0  # 1/4弦长位置
        self.dY = -1.0  # Y方向偏移量
        self.dZ = 0  # Z方向偏移量

        # derived description
        self.S = 0  # 参考机翼面积
        self.AR = 0  # 展弦比
        self.A_25 = 0  # 1/4弦线后掠角
        self.TaperRatio = 0  # 梯形比
        self.MAC = 0  # 平均气动弦长

        # widgets
        self.widget = QWidget()
        self.airfoil_label = QLabel('翼型:')
        self.airfoil_combobox = QComboBox()
        self.span_label = QLabel('展长(m):')
        self.span_lineedit = QLineEdit()
        self.section_num_label = QLabel('控制剖面数量:')
        self.section_num_lineedit = QLineEdit()
        self.wing_root_len_label = QLabel('翼根长度(m):')
        self.wing_root_len_lineedit = QLineEdit()
        self.wing_tip_len_label = QLabel('翼尖长度(m):')
        self.wing_tip_len_lineedit = QLineEdit()
        self.sweep_back_label = QLabel('前缘后掠角(°):')
        self.sweep_back_lineedit = QLineEdit()
        self.dihedral_label = QLabel('上反角(°):')
        self.dihedral_lineedit = QLineEdit()
        self.twist_label = QLabel('扭转角(°):')
        self.twist_lineedit = QLineEdit()
        self.x25_label = QLabel('1/4弦线X轴位置(m):')
        self.x25_lineedit = QLineEdit()
        self.dy_label = QLabel('Y方向偏移量(m):')
        self.dy_lineedit = QLineEdit()
        self.dz_label = QLabel('Z方向偏移量(m):')
        self.dz_lineedit = QLineEdit()

        self.ref_area = QLabel('参考面积: %.2f' % self.S)
        self.aspect_ratio = QLabel('展弦比: %.2f' % self.AR)
        self.sweep_back_25 = QLabel('1/4弦线后掠角: %.2f' % self.A_25)
        self.taper_ratio = QLabel('梯形比: %.2f' % self.TaperRatio)
        self.mac = QLabel('平均气动弦长: %.2f' % self.MAC)

        # layout
        self.layout = QVBoxLayout()
        self.label = QLabel('机翼设计参数:')
        self.installation_label = QLabel('安装参数:')
        self.param_layout = QHBoxLayout()
        self.intrinsic_param_layout = QGridLayout()
        self.derived_param_layout = QVBoxLayout()

    def init_widget(self):
        self.widget.setLayout(self.layout)
        self.layout.addWidget(self.label)
        self.layout.addLayout(self.param_layout)
        self.param_layout.addLayout(self.intrinsic_param_layout)
        self.param_layout.addLayout(self.derived_param_layout)

        self.airfoil_combobox.addItems(Wing.airfoil_list)
        self.intrinsic_param_layout.addWidget(self.airfoil_label, 0, 0)
        self.intrinsic_param_layout.addWidget(self.airfoil_combobox, 0, 1)

        self.span_lineedit.setText(str(self.Span))
        self.intrinsic_param_layout.addWidget(self.span_label, 1, 0)
        self.intrinsic_param_layout.addWidget(self.span_lineedit, 1, 1)

        self.section_num_lineedit.setText(str(self.SectionNum))
        self.intrinsic_param_layout.addWidget(self.section_num_label, 2, 0)
        self.intrinsic_param_layout.addWidget(self.section_num_lineedit, 2, 1)

        self.wing_root_len_lineedit.setText(str(self.C_root))
        self.intrinsic_param_layout.addWidget(self.wing_root_len_label, 3, 0)
        self.intrinsic_param_layout.addWidget(self.wing_root_len_lineedit, 3, 1)

        self.wing_tip_len_lineedit.setText(str(self.C_tip))
        self.intrinsic_param_layout.addWidget(self.wing_tip_len_label, 4, 0)
        self.intrinsic_param_layout.addWidget(self.wing_tip_len_lineedit, 4, 1)

        self.sweep_back_lineedit.setText(str(self.SweepBack))
        self.intrinsic_param_layout.addWidget(self.sweep_back_label, 5, 0)
        self.intrinsic_param_layout.addWidget(self.sweep_back_lineedit, 5, 1)

        self.dihedral_lineedit.setText(str(self.Dihedral))
        self.intrinsic_param_layout.addWidget(self.dihedral_label, 6, 0)
        self.intrinsic_param_layout.addWidget(self.dihedral_lineedit, 6, 1)

        self.twist_lineedit.setText(str(self.Twist))
        self.intrinsic_param_layout.addWidget(self.twist_label, 7, 0)
        self.intrinsic_param_layout.addWidget(self.twist_lineedit, 7, 1)

        self.intrinsic_param_layout.addWidget(self.installation_label, 8, 0)

        self.x25_lineedit.setText(str(self.X_25))
        self.intrinsic_param_layout.addWidget(self.x25_label, 9, 0)
        self.intrinsic_param_layout.addWidget(self.x25_lineedit, 9, 1)

        self.dy_lineedit.setText(str(self.dY))
        self.intrinsic_param_layout.addWidget(self.dy_label, 10, 0)
        self.intrinsic_param_layout.addWidget(self.dy_lineedit, 10, 1)

        self.dz_lineedit.setText(str(self.dZ))
        self.intrinsic_param_layout.addWidget(self.dz_label, 11, 0)
        self.intrinsic_param_layout.addWidget(self.dz_lineedit, 11, 1)

        self.derived_param_layout.addWidget(self.ref_area)
        self.derived_param_layout.addWidget(self.aspect_ratio)
        self.derived_param_layout.addWidget(self.sweep_back_25)
        self.derived_param_layout.addWidget(self.taper_ratio)
        self.derived_param_layout.addWidget(self.mac)

    def update_derived_param(self):
        # get intrinsic param
        self.Airfoil = Wing.airfoil_dir + self.airfoil_combobox.currentText()
        self.Span = float(self.span_lineedit.text())
        self.SectionNum = int(self.section_num_lineedit.text())
        self.C_root = float(self.wing_root_len_lineedit.text())
        self.C_tip = float(self.wing_tip_len_lineedit.text())
        self.SweepBack = float(self.sweep_back_lineedit.text())
        self.Dihedral = float(self.dihedral_lineedit.text())
        self.Twist = float(self.twist_lineedit.text())

        self.X_25 = float(self.x25_lineedit.text())
        self.dY = float(self.dy_lineedit.text())
        self.dZ = float(self.dz_lineedit.text())

        # calc derived param
        self.S = float(0.5 * self.Span * (self.C_root + self.C_tip))
        self.AR = float(math.pow(self.Span, 2) / self.S)
        self.TaperRatio = float(self.C_tip / self.C_root)
        self.A_25 = math.degrees(math.atan(
            math.tan(math.radians(self.SweepBack)) - (1 - self.TaperRatio) / (self.TaperRatio * (1 + self.TaperRatio))))
        self.MAC = float(2 / 3 * self.C_root * (1 - math.pow(self.TaperRatio, 3)) / (1 - math.pow(self.TaperRatio, 2)))

        self.ref_area.setText('参考面积: %.2f' % self.S)
        self.aspect_ratio.setText('展弦比: %.2f' % self.AR)
        self.sweep_back_25.setText('1/4弦线后掠角: %.2f' % self.A_25)
        self.taper_ratio.setText('梯形比: %.2f' % self.TaperRatio)
        self.mac.setText('平均气动弦长: %.2f' % self.MAC)

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
        self.X_25 = 52.0
        self.dY = 0
        self.dZ = 0

        # derived description
        self.V_v = 0.079  # 垂尾容量

        # widget
        self.label.setText('垂尾设计参数:')
        self.capacity = QLabel('垂尾容量: %.4f' % self.V_v)

        # layout
        self.derived_param_layout.addWidget(self.capacity)

    def update_derived_param(self):
        Wing.update_derived_param(self)
        self.capacity.setText('垂尾容量: %.4f' % self.V_v)


# 平尾
class HorizontalStabilizer(Wing):
    def __init__(self, _filename):
        Wing.__init__(self, _filename)

        # real description
        self.SectionNum = 4

        # intrinsic description
        self.Span = 10.2
        self.C_root = 5.0
        self.C_tip = 1.9
        self.SweepBack = 25
        self.Dihedral = 0
        self.Twist = 0

        # installation description
        self.X_25 = 54.0
        self.dY = 0
        self.dZ = 0

        # derived description
        self.V_h = 0

        # widget
        self.label.setText('平尾设计参数:')
        self.capacity = QLabel('平尾容量: %.4f' % self.V_h)

        # layout
        self.derived_param_layout.addWidget(self.capacity)

    def update_derived_param(self):
        Wing.update_derived_param(self)
        self.capacity.setText('平尾容量: %.4f' % self.V_h)
