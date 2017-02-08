import matplotlib.pyplot as plt
import numpy as np
import math
from stl import mesh
from PyQt5.QtWidgets import *
import os


class Airfoil(object):
    '''2维翼型描述, 弦长为单位1'''

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


class Section(object):
    ''''组成3维机翼的真实剖面'''

    def __init__(self, _airfoil, _z, _x_front, _x_tail, _dy, _theta):
        '''
        :param _airfoil: 剖面上的翼型
        :param _z: 剖面与根部的距离
        :param _x_front: 前缘点在x方向上的位置
        :param _x_tail: 后缘点在x方向上的位置
        :param _dy: 由于机翼上反角导致的在y方向上引起的偏移
        :param _theta: 绕前缘的扭转角
        '''

        self.airfoil = _airfoil
        self.z = _z
        self.x_front = _x_front
        self.x_tail = _x_tail
        self.dy = _dy
        self.theta = _theta

        self.n = len(_airfoil.x)
        self.chord = np.zeros((self.n, 3))  # 剖面弦线上的离散点
        self.up_pts = np.zeros((self.n, 3))  # 剖面上表面的点
        self.down_pts = np.zeros((self.n, 3))  # 剖面下表面的点

    def calc_discrete_pts(self):
        '''计算描述剖面形状的空间坐标点'''

        # 弦长
        chord_len = self.x_tail - self.x_front

        # 计算正常放置时上、下表面上每点的坐标
        for i in range(0, self.n):
            '''x方向'''
            self.up_pts[i][0] = self.x_front + chord_len * self.airfoil.x[i]
            self.down_pts[i][0] = self.up_pts[i][0]
            '''y方向'''
            self.up_pts[i][1] = chord_len * self.airfoil.y_up[i]
            self.down_pts[i][1] = chord_len * self.airfoil.y_down[i]
            '''z方向'''
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


class Wing(object):
    '''简单机翼，由一个个剖面组成'''

    '翼型列表'
    airfoil_list = []
    airfoil_dir = './airfoil/'

    @staticmethod
    def update_airfoil_list():
        for f in os.listdir(Wing.airfoil_dir):
            cur_filename = os.path.join(Wing.airfoil_dir, f)
            if os.path.isfile(cur_filename) and cur_filename.endswith('.dat'):
                Wing.airfoil_list.append(f)

    def __init__(self, _filename):
        self.filename = _filename
        self.SectionNum = 12
        self.section_list = []
        self.widget = QWidget()

        'intrinsic description'
        self.Airfoil = ''  # 翼型
        self.Span = 6380  # 展长
        self.C_root = 620  # 翼根弦长
        self.C_tip = 258  # 翼尖弦长
        self.SweepBack = 30  # 后掠角
        self.Dihedral = 5  # 上反角
        self.Twist = -4  # 绕前缘扭转角
        self.X_25 = 1600  # 1/4弦长位置

        'derived description'
        self.S = 0  # 参考机翼面积
        self.AR = 0  # 展弦比
        self.A_25 = 0  # 1/4弦线后掠角
        self.TaperRatio = 0  # 梯形比
        self.MAC = 0  # 平均气动弦长

    def init_widget(self, _widget_name):
        layout = QVBoxLayout()
        self.widget.setLayout(layout)

        label = QLabel('机翼设计参数：')
        layout.addWidget(label)

        param_layout = QHBoxLayout()
        layout.addLayout(param_layout)

        intrinsic_param_layout = QGridLayout()
        derived_param_layout = QVBoxLayout()

        param_layout.addLayout(intrinsic_param_layout)
        param_layout.addLayout(derived_param_layout)

        airfoil_label = QLabel('翼型:')
        self.airfoil_combobox = QComboBox()
        self.airfoil_combobox.addItems(Wing.airfoil_list)
        intrinsic_param_layout.addWidget(airfoil_label, 0, 0)
        intrinsic_param_layout.addWidget(self.airfoil_combobox, 0, 1)

        span_label = QLabel('展长(cm):')
        self.span_lineedit = QLineEdit()
        self.span_lineedit.setText(str(self.Span))
        intrinsic_param_layout.addWidget(span_label, 1, 0)
        intrinsic_param_layout.addWidget(self.span_lineedit, 1, 1)

        section_num_label = QLabel('控制剖面数量:')
        self.section_num_lineedit = QLineEdit()
        self.section_num_lineedit.setText(str(self.SectionNum))
        intrinsic_param_layout.addWidget(section_num_label, 2, 0)
        intrinsic_param_layout.addWidget(self.section_num_lineedit, 2, 1)

        wing_root_len_label = QLabel('翼根长度(cm):')
        self.wing_root_len_lineedit = QLineEdit()
        self.wing_root_len_lineedit.setText(str(self.C_root))
        intrinsic_param_layout.addWidget(wing_root_len_label, 3, 0)
        intrinsic_param_layout.addWidget(self.wing_root_len_lineedit, 3, 1)

        wing_tip_len_label = QLabel('翼尖长度(cm):')
        self.wing_tip_len_lineedit = QLineEdit()
        self.wing_tip_len_lineedit.setText(str(self.C_tip))
        intrinsic_param_layout.addWidget(wing_tip_len_label, 4, 0)
        intrinsic_param_layout.addWidget(self.wing_tip_len_lineedit, 4, 1)

        sweep_back_label = QLabel('前缘后掠角(°):')
        self.sweep_back_lineedit = QLineEdit()
        self.sweep_back_lineedit.setText(str(self.SweepBack))
        intrinsic_param_layout.addWidget(sweep_back_label, 5, 0)
        intrinsic_param_layout.addWidget(self.sweep_back_lineedit, 5, 1)

        dihedral_label = QLabel('上反角(°):')
        self.dihedral_lineedit = QLineEdit()
        self.dihedral_lineedit.setText(str(self.Dihedral))
        intrinsic_param_layout.addWidget(dihedral_label, 6, 0)
        intrinsic_param_layout.addWidget(self.dihedral_lineedit, 6, 1)

        twist_label = QLabel('扭转角(°):')
        self.twist_lineedit = QLineEdit()
        self.twist_lineedit.setText(str(self.Twist))
        intrinsic_param_layout.addWidget(twist_label, 7, 0)
        intrinsic_param_layout.addWidget(self.twist_lineedit, 7, 1)

        x25_label = QLabel('1/4弦线X轴位置(cm):')
        self.x25_lineedit = QLineEdit()
        self.x25_lineedit.setText(str(self.X_25))
        intrinsic_param_layout.addWidget(x25_label, 8, 0)
        intrinsic_param_layout.addWidget(self.x25_lineedit, 8, 1)

        self.ref_area = QLabel('参考面积: %.2f' % (self.S))
        self.aspect_ratio = QLabel('展弦比: %.2f' % (self.AR))
        self.sweep_back_25 = QLabel('1/4弦线后掠角: %.2f' % (self.A_25))
        self.taper_ratio = QLabel('梯形比: %.2f' % (self.TaperRatio))
        self.mac = QLabel('平均气动弦长: %.2f' % (self.MAC))
        derived_param_layout.addWidget(self.ref_area)
        derived_param_layout.addWidget(self.aspect_ratio)
        derived_param_layout.addWidget(self.sweep_back_25)
        derived_param_layout.addWidget(self.taper_ratio)
        derived_param_layout.addWidget(self.mac)

        return self.widget, _widget_name

    def get_intrinsic_param(self):
        self.Airfoil = Wing.airfoil_dir + self.airfoil_combobox.currentText()
        self.Span = float(self.span_lineedit.text())
        self.SectionNum = int(self.section_num_lineedit.text())
        self.C_root = float(self.wing_root_len_lineedit.text())
        self.C_tip = float(self.wing_tip_len_lineedit.text())
        self.SweepBack = float(self.sweep_back_lineedit.text())
        self.Dihedral = float(self.dihedral_lineedit.text())
        self.Twist = float(self.twist_lineedit.text())
        self.X_25 = float(self.x25_lineedit.text())

    def update_derived_param(self):
        self.get_intrinsic_param()

        self.S = float(0.5 * self.Span * (self.C_root + self.C_tip))
        self.AR = float(math.pow(self.Span, 2) / self.S)
        self.TaperRatio = float(self.C_tip / self.C_root)
        self.A_25 = math.degrees(math.atan(math.tan(math.radians(self.SweepBack)) - (1 - self.TaperRatio) / (
            self.TaperRatio * (1 + self.TaperRatio))))
        self.MAC = float(
            2 / 3 * self.C_root * (1 - math.pow(self.TaperRatio, 3)) / (1 - math.pow(self.TaperRatio, 2)))

        self.ref_area.setText('参考面积: %.2f' % (self.S))
        self.aspect_ratio.setText('展弦比: %.2f' % (self.AR))
        self.sweep_back_25.setText('1/4弦线后掠角: %.2f' % (self.A_25))
        self.taper_ratio.setText('梯形比: %.2f' % (self.TaperRatio))
        self.mac.setText('平均气动弦长: %.2f' % (self.MAC))

    def set_real_description(self, _airfoil, _z, _x_front, _x_tail, _dy, _theta):
        self.airfoil_list = _airfoil
        self.z_list = _z
        self.x_front_list = _x_front
        self.x_tail_list = _x_tail
        self.dy_list = _dy
        self.theta_list = _theta

        self.SectionNum = len(_z)

    def calc_sections(self):
        for i in range(0, self.SectionNum):
            self.section_list.append(
                Section(self.airfoil_list[i], self.z_list[i], self.x_front_list[i], self.x_tail_list[i],
                        self.dy_list[i], self.theta_list[i]))
        for i in range(0, len(self.z_list)):
            self.section_list[i].calc_discrete_pts()

    def generate_wing_stl(self):
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


class VerticalStabilizer(Wing):
   '''垂尾'''

   def __init__(self):
       pass