import sys
import os
import math
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from src.aircraft.wing import *

airfoil_dir = '../../airfoil/'  # 翼型文件夹路径
airfoil_list = []  # 翼型列表
z = np.array([0., 0.624029, 1.38967, 2.43503, 3.73439, 5.25574, 6.96162, 8.81003, 10.7555, 12.75, 14.7445, 16.69, 18.5384, 20.2443, 21.7656, 23.065, 24.1103, 24.876, 25.343, 25.5])
x_front = np.array([0., 0.05, 0.3, 1.7, 4.9, 6.85, 8.45, 9.65, 10.6, 11.1, 11.7, 12.1, 12.4, 12.8, 13.2, 13.7, 14.1, 14.5, 15.2, 16.])
y_front = np.zeros(len(z), float)
x_tail = np.array([19.7, 19.6, 19.6, 19.5, 19.3, 19, 18.3, 17.3, 16.6, 16.5, 16.8, 17, 17.45, 17.8, 18.1, 18.4, 18.55, 18.65, 18.3, 17.8])
y_tail = np.zeros(len(z), float)
thickness = np.ones(len(z), float)
airfoil = []

for i in range(0, len(z)):
    airfoil.append(Airfoil("../../airfoil/naca0012.dat"))


def update_airfoil_list():
    for f in os.listdir(airfoil_dir):
        cur_filename = os.path.join(airfoil_dir, f)
        if os.path.isfile(cur_filename) and cur_filename.endswith('.dat'):
            airfoil_list.append(f)


class BWB_Section(object):
    # 描述机翼的剖面
    BWB_SEC_PARAM = ['Airfoil', 'Thickness Ratio', 'Z(m)', 'X_front(m)', 'Y_front(m)', 'X_tail(m)', 'Y_tail(m)']

    def __init__(self, airfoil, thickness, z, xf, yf, xt, yt):
        # Descriptions
        self.airfoil = airfoil
        self.thickness = thickness
        self.z = z
        self.xf = xf
        self.yf = yf
        self.xt = xt
        self.yt = yt

        # Widgets
        self.airfoil_combobox = QComboBox()
        self.airfoil_combobox.addItems(airfoil_list)
        self.thick_le = QLineEdit()
        self.thick_le.setText(str(self.thickness))
        self.z_le = QLineEdit()
        self.z_le.setText(str(self.z))
        self.xf_le = QLineEdit()
        self.xf_le.setText(str(self.xf))
        self.yf_le = QLineEdit()
        self.yf_le.setText(str(self.yf))
        self.xt_le = QLineEdit()
        self.xt_le.setText(str(self.xt))
        self.yt_le = QLineEdit()
        self.yt_le.setText(str(self.yt))

        self.widget_list = []
        self.widget_list.append(self.airfoil_combobox)
        self.widget_list.append(self.thick_le)
        self.widget_list.append(self.z_le)
        self.widget_list.append(self.xf_le)
        self.widget_list.append(self.yf_le)
        self.widget_list.append(self.xt_le)
        self.widget_list.append(self.yt_le)


class BWB_UI(QDialog):
    def __init__(self, _airfoil_list, _tck_list, _z_list, _xf_list, _yf_list, _xt_list, _yt_list):
        super(BWB_UI, self).__init__()

        if len(_airfoil_list) != len(_z_list):
            raise Exception("Invalid Parameter!")
        if len(_z_list) != len(_xf_list):
            raise Exception("Invalid Parameter!")
        if len(_xf_list) != len(_yf_list):
            raise Exception("Invalid Parameter!")
        if len(_yf_list) != len(_xt_list):
            raise Exception("Invalid Parameter!")
        if len(_xt_list) != len(_yt_list):
            raise Exception("Invalid Parameter!")
        if len(_yt_list) != len(_tck_list):
            raise Exception("Invalid Parameter!")

        # Contents
        self.row_cnt = len(_airfoil_list)
        self.col_cnt = len(BWB_Section.BWB_SEC_PARAM)
        self.section = []

        self.airfoil_list = _airfoil_list
        self.thickness_list = _tck_list
        self.z_list = _z_list
        self.xf_list = _xf_list
        self.yf_list = _yf_list
        self.xt_list = _xt_list
        self.yt_list = _yt_list

        # Param Table
        self.section_param_list = QTableWidget(self.row_cnt, self.col_cnt)
        self.section_param_list.setHorizontalHeaderLabels(BWB_Section.BWB_SEC_PARAM)
        self.section_param_list.horizontalHeader().setStretchLastSection(True)  # 让表格填满窗口
        self.section_param_list.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        # Buttons
        self.btn_model = QPushButton()
        self.btn_model.setText('Get IGES Model')
        self.btn_model.clicked.connect(self.gen_model)

        self.btn_mesh = QPushButton()
        self.btn_mesh.setText('Get Plot3D Mesh')
        self.btn_mesh.clicked.connect(self.gen_mesh)

        self.btn_exit = QPushButton()
        self.btn_exit.setText('Quit')
        self.btn_exit.clicked.connect(self.close)

        # Layout
        self.global_layout = QVBoxLayout()
        self.param_layout = QHBoxLayout()
        self.btn_layout = QHBoxLayout()
        self.param_layout.addWidget(self.section_param_list)
        self.btn_layout.addWidget(self.btn_model)
        self.btn_layout.addWidget(self.btn_mesh)
        self.btn_layout.addWidget(self.btn_exit)

        # Base Frame
        self.setWindowTitle('BWB Parametric Modeling')
        self.setWindowIcon(QIcon('../../pic/plane_icon.jpg'))
        self.setWindowFlags(Qt.WindowMaximizeButtonHint | Qt.WindowMinimizeButtonHint | Qt.WindowCloseButtonHint)
        self.resize(750, 450)
        self.setLayout(self.global_layout)
        self.global_layout.addLayout(self.param_layout)
        self.global_layout.addLayout(self.btn_layout)

        # Init Contents
        for i in range(0, self.row_cnt):
            csec = BWB_Section(self.airfoil_list[i], self.thickness_list[i], self.z_list[i], self.xf_list[i], self.yf_list[i], self.xt_list[i], self.yt_list[i])
            self.section.append(csec)
            for j in range(0, self.col_cnt):
                self.section_param_list.setCellWidget(i, j, csec.widget_list[j])

    def update_content(self):
        pass

    def gen_model(self):
        pass

    def gen_mesh(self):
        pass


if __name__ == '__main__':
    update_airfoil_list()
    app = QApplication(sys.argv)
    pm = BWB_UI(airfoil, thickness, z, x_front, y_front, x_tail, y_tail)
    pm.show()
    sys.exit(app.exec_())
