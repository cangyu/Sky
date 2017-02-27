from PyQt5.QtWidgets import *
from aeroplane import Aircraft
import sys
import os
import math

airfoil_dir = './airfoil/'
airfoil_list = []  # 翼型列表


def update_airfoil_list():
    for f in os.listdir(airfoil_dir):
        cur_filename = os.path.join(airfoil_dir, f)
        if os.path.isfile(cur_filename) and cur_filename.endswith('.dat'):
            airfoil_list.append(f)


class WingUI(QWidget):
    def __init__(self):
        QWidget.__init__(self)

        # intrinsic description
        self.Airfoil = ''  # 翼型
        self.Span = 30  # 展长
        self.SectionNum = 12  # 剖面数量
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

        # widgets
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
        self.global_layout = QHBoxLayout()
        self.intrinsic_param_layout = QGridLayout()
        self.ipl_cnt = 0

        self.derived_param_layout = QGridLayout()
        self.dpl_cnt = 0

    def initUI(self):
        # global layout configuration
        self.setLayout(self.global_layout)

        self.global_layout.addLayout(self.intrinsic_param_layout)
        self.global_layout.addLayout(self.derived_param_layout)

        # intrinsic params part
        self.airfoil_combobox.addItems(airfoil_list)
        self.intrinsic_param_layout.addWidget(self.airfoil_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.airfoil_combobox, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.span_lineedit.setText(str(self.Span))
        self.intrinsic_param_layout.addWidget(self.span_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.span_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.section_num_lineedit.setText(str(self.SectionNum))
        self.intrinsic_param_layout.addWidget(self.section_num_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.section_num_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.wing_root_len_lineedit.setText(str(self.C_root))
        self.intrinsic_param_layout.addWidget(self.wing_root_len_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.wing_root_len_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.wing_tip_len_lineedit.setText(str(self.C_tip))
        self.intrinsic_param_layout.addWidget(self.wing_tip_len_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.wing_tip_len_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.sweep_back_lineedit.setText(str(self.SweepBack))
        self.intrinsic_param_layout.addWidget(self.sweep_back_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.sweep_back_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.dihedral_lineedit.setText(str(self.Dihedral))
        self.intrinsic_param_layout.addWidget(self.dihedral_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.dihedral_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.twist_lineedit.setText(str(self.Twist))
        self.intrinsic_param_layout.addWidget(self.twist_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.twist_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.x25_lineedit.setText(str(self.X_25))
        self.intrinsic_param_layout.addWidget(self.x25_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.x25_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.dy_lineedit.setText(str(self.dY))
        self.intrinsic_param_layout.addWidget(self.dy_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.dy_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.dz_lineedit.setText(str(self.dZ))
        self.intrinsic_param_layout.addWidget(self.dz_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.dz_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        # derived params part
        self.derived_param_layout.addWidget(self.ref_area, self.dpl_cnt, 0)
        self.dpl_cnt += 1

        self.derived_param_layout.addWidget(self.aspect_ratio, self.dpl_cnt, 0)
        self.dpl_cnt += 1

        self.derived_param_layout.addWidget(self.sweep_back_25, self.dpl_cnt, 0)
        self.dpl_cnt += 1

        self.derived_param_layout.addWidget(self.taper_ratio, self.dpl_cnt, 0)
        self.dpl_cnt += 1

        self.derived_param_layout.addWidget(self.mac, self.dpl_cnt, 0)
        self.dpl_cnt += 1

    def updateUI(self):
        # get intrinsic param
        self.Airfoil = airfoil_dir + self.airfoil_combobox.currentText()
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


class VerticalStabilizerUI(WingUI):
    def __init__(self):
        WingUI.__init__(self)

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
        self.V_v = 0  # 垂尾容量

        # widget
        self.capacity_label = QLabel('垂尾容量: %.4f' % self.V_v)

    def initUI(self):
        WingUI.initUI(self)

        self.derived_param_layout.addWidget(self.capacity_label, self.dpl_cnt, 0)
        self.dpl_cnt += 1

    def updateUI(self):
        WingUI.initUI(self)

        self.capacity_label.setText('垂尾容量: %.4f' % self.V_v)


class HorizontalStabilizerUI(WingUI):
    def __init__(self):
        WingUI.__init__(self)

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

        # widget
        self.capacity_label = QLabel('平尾容量: %.4f' % self.V_h)

    def initUI(self):
        WingUI.initUI(self)

        self.derived_param_layout.addWidget(self.capacity_label, self.dpl_cnt, 0)
        self.dpl_cnt += 1

    def updateUI(self):
        WingUI.initUI(self)

        self.capacity_label.setText('平尾容量: %.4f' % self.V_h)


class FuselageUI(QWidget):
    def __init__(self):
        QWidget.__init__(self)

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

        # widget
        self.L_label = QLabel('机身全长(m):')
        self.L_lineedit = QLineEdit()
        self.L_lineedit.setText(str(self.L))

        self.D_label = QLabel('机身直径(m):')
        self.D_lineedit = QLineEdit()
        self.D_lineedit.setText(str(self.D))

        self.Theta_fc_label = QLabel('擦地角(°):')
        self.Theta_fc_lineedit = QLineEdit()
        self.Theta_fc_lineedit.setText(str(self.Theta_fc))

        self.r0_label = QLabel('机头占全长比例:')
        self.r0_dsb = QDoubleSpinBox()
        self.r0_dsb.setRange(0.0, 1.0)
        self.r0_dsb.setSingleStep(0.01)
        self.r0_dsb.setValue(self.r0)

        self.r2_label = QLabel('机尾占全长比例:')
        self.r2_dsb = QDoubleSpinBox()
        self.r2_dsb.setRange(0.0, 1.0)
        self.r2_dsb.setSingleStep(0.01)
        self.r2_dsb.setValue(self.r2)

        self.LvD_total_label = QLabel('机身长径比: %.2f' % self.lambda_total)
        self.LvD_front_label = QLabel('头部长径比: %.2f' % self.lambda_front)
        self.LvD_mid_label = QLabel('中部长径比: %.2f' % self.lambda_mid)
        self.LvD_tail_label = QLabel('尾部长径比: %.2f' % self.lambda_tail)
        self.HeadingDirCapacity_label = QLabel('航向机身容量: %.2f' % self.HeadingDirCapacity)
        self.PitchDirCapacity_label = QLabel('纵向机身容量: %.2f' % self.PitchDirCapacity)

        # layout
        self.global_layout = QHBoxLayout()

        self.intrinsic_param_layout = QGridLayout()
        self.ipl_cnt = 0

        self.derived_param_layout = QGridLayout()
        self.dpl_cnt = 0

    def initUI(self):
        # global layout
        self.setLayout(self.global_layout)
        self.global_layout.addLayout(self.intrinsic_param_layout)
        self.global_layout.addLayout(self.derived_param_layout)

        # intrinsic params part
        self.intrinsic_param_layout.addWidget(self.L_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.L_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.intrinsic_param_layout.addWidget(self.D_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.D_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.intrinsic_param_layout.addWidget(self.Theta_fc_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.Theta_fc_lineedit, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.intrinsic_param_layout.addWidget(self.r0_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.r0_dsb, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        self.intrinsic_param_layout.addWidget(self.r2_label, self.ipl_cnt, 0)
        self.intrinsic_param_layout.addWidget(self.r2_dsb, self.ipl_cnt, 1)
        self.ipl_cnt += 1

        # derived params part
        self.derived_param_layout.addWidget(self.LvD_total_label, self.dpl_cnt, 0)
        self.dpl_cnt += 1

        self.derived_param_layout.addWidget(self.LvD_front_label, self.dpl_cnt, 0)
        self.dpl_cnt += 1

        self.derived_param_layout.addWidget(self.LvD_mid_label, self.dpl_cnt, 0)
        self.dpl_cnt += 1

        self.derived_param_layout.addWidget(self.LvD_tail_label, self.dpl_cnt, 0)
        self.dpl_cnt += 1

        self.derived_param_layout.addWidget(self.HeadingDirCapacity_label, self.dpl_cnt, 0)
        self.dpl_cnt += 1

        self.derived_param_layout.addWidget(self.PitchDirCapacity_label, self.dpl_cnt, 0)
        self.dpl_cnt += 1

    def updateUI(self):
        # get new intrinsic param
        self.L = float(self.L_lineedit.text())
        self.D = float(self.D_lineedit.text())
        self.Theta_fc = float(self.Theta_fc_lineedit.text())
        self.r0 = self.r0_dsb.value()
        self.r2 = self.r2_dsb.value()

        # calc derived param
        self.r1 = 1.0 - self.r0 - self.r2
        self.L0 = self.L * self.r0
        self.L1 = self.L * self.r1
        self.L2 = self.L * self.r2

        self.lambda_total = self.L / self.D
        self.lambda_front = self.L0 / self.D
        self.lambda_mid = self.L1 / self.D
        self.lambda_tail = self.L2 / self.D

        # show changed
        self.LvD_total_label.setText('机身长径比: %.6f' % self.lambda_total)
        self.LvD_front_label.setText('头部长径比: %.6f' % self.lambda_front)
        self.LvD_mid_label.setText('中部长径比: %.6f' % self.lambda_mid)
        self.LvD_tail_label.setText('尾部长径比: %.6f' % self.lambda_tail)
        self.HeadingDirCapacity_label.setText('航向机身容量: %.6f' % self.HeadingDirCapacity)
        self.PitchDirCapacity_label.setText('纵向机身容量: %.6f' % self.PitchDirCapacity)


class AircraftUI(QWidget):
    def __init__(self):
        QWidget.__init__(self)

        # entity
        self.aircraft = Aircraft()

        # components
        self.wingUI = WingUI()
        self.hsUI = HorizontalStabilizerUI()
        self.vsUI = VerticalStabilizerUI()
        self.fuselageUI = FuselageUI()

        # widget
        self.comp_tab = QTabWidget()

        self.btn_calc = QPushButton()
        self.btn_calc.setText('计算')
        self.btn_calc.setToolTip('根据原始参数计算衍生参数')
        self.btn_calc.clicked.connect(self.update_param)

        self.btn_gen = QPushButton()
        self.btn_gen.setText('生成')
        self.btn_gen.setToolTip('根据输入参数生成模型')
        self.btn_gen.clicked.connect(self.gen_model)

        self.btn_exit = QPushButton()
        self.btn_exit.setText('退出')
        self.btn_exit.setToolTip('退出本程序')
        self.btn_exit.clicked.connect(self.close)

        # layout
        self.global_layout = QVBoxLayout()

        self.param_layout = QVBoxLayout()
        self.btn_layout = QHBoxLayout()

        update_airfoil_list()

    def initUI(self):
        self.setWindowTitle('简单客机一体化设计')

        # global layout
        self.setLayout(self.global_layout)
        self.global_layout.addLayout(self.param_layout)
        self.global_layout.addLayout(self.btn_layout)

        # params part
        self.param_layout.addWidget(self.comp_tab)

        self.comp_tab.addTab(self.fuselageUI, '机身')
        self.comp_tab.addTab(self.wingUI, '机翼')
        self.comp_tab.addTab(self.vsUI, '垂尾')
        self.comp_tab.addTab(self.hsUI, '平尾')

        # button part
        self.btn_layout.addWidget(self.btn_calc)
        self.btn_layout.addWidget(self.btn_gen)
        self.btn_layout.addWidget(self.btn_exit)

        # initialize component ui
        self.wingUI.initUI()
        self.fuselageUI.initUI()
        self.hsUI.initUI()
        self.vsUI.initUI()

    def update_param(self):
        self.aircraft.update_derived_param(self)

    def gen_model(self):
        self.aircraft.generate(self)


if __name__ == '__main__':
    app = QApplication(sys.argv)

    ythsj = AircraftUI()
    ythsj.initUI()
    ythsj.show()

    sys.exit(app.exec_())
