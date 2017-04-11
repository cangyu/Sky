from PyQt5.QtWidgets import *
from PyQt5.QtGui import QIcon
import sys
import os
import math
from src.aircraft.wing import *

# 翼型文件夹路径
airfoil_dir = '../airfoil/'

# 翼型列表
airfoil_list = []

z = np.array([0., 0.624029, 1.38967, 2.43503, 3.73439, 5.25574, 6.96162, 8.81003, 10.7555, 12.75, 14.7445, 16.69, 18.5384, 20.2443, 21.7656, 23.065, 24.1103, 24.876, 25.343, 25.5])
x_front = np.array([0., 0.05, 0.3, 1.7, 4.9, 6.85, 8.45, 9.65, 10.6, 11.1, 11.7, 12.1, 12.4, 12.8, 13.2, 13.7, 14.1, 14.5, 15.2, 16.])
y_front = np.linspace(0, 0, len(z))
x_tail = np.array([19.7, 19.6, 19.6, 19.5, 19.3, 19, 18.3, 17.3, 16.6, 16.5, 16.8, 17, 17.45, 17.8, 18.1, 18.4, 18.55, 18.65, 18.3, 17.8])
y_tail = np.linspace(0, 0, len(z))
thickness = np.ones(len(z))
airfoil = []
for i in range(0, len(z)):
    airfoil.append(Airfoil("../airfoil/naca0012.dat"))


def update_airfoil_list():
    for f in os.listdir(airfoil_dir):
        cur_filename = os.path.join(airfoil_dir, f)
        if os.path.isfile(cur_filename) and cur_filename.endswith('.dat'):
            airfoil_list.append(f)


class Section(QWidget):
    def __init__(self, airfoil, z, xf, yf, xt, yt, thickness):
        super(Section, self).__init__()

        # Descriptions
        self.airfoil = airfoil
        self.z = z
        self.xf = xf
        self.yf = yf
        self.xt = xt
        self.yt = yt
        self.thickness = thickness

        # Layout
        self.global_layout = QHBoxLayout()

        # Widgets
        self.airfoil_combobox = QComboBox()
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
        self.thick_le = QLineEdit()
        self.thick_le.setText(str(self.thickness))

        self.airfoil_combobox.addItems(airfoil_list)
        self.setLayout(self.global_layout)

        # Init
        self.global_layout.addWidget(self.airfoil_combobox)
        self.global_layout.addWidget(self.z_le)
        self.global_layout.addWidget(self.xf_le)
        self.global_layout.addWidget(self.yf_le)
        self.global_layout.addWidget(self.xt_le)
        self.global_layout.addWidget(self.yt_le)
        self.global_layout.addWidget(self.thick_le)


class BWB_UI(QWidget):
    def __init__(self, _airfoil_list, _z_list, _xf_list, _yf_list, _xt_list, _yt_list, _tck_list):
        QWidget.__init__(self)

        assert len(_airfoil_list) == len(_z_list)
        assert len(_z_list) == len(_xf_list)
        assert len(_xf_list) == len(_yf_list)
        assert len(_yf_list) == len(_xt_list)
        assert len(_xt_list) == len(_yt_list)
        assert len(_yt_list) == len(_tck_list)

        # Contents
        self.sec_cnt = len(_airfoil_list)
        self.section = []
        self.airfoil_list = _airfoil_list
        self.z_list = _z_list
        self.xf_list = _xf_list
        self.yf_list = _yf_list
        self.xt_list = _xt_list
        self.yt_list = _yt_list
        self.thickness_list = _tck_list

        # Layout
        self.global_layout = QVBoxLayout()
        self.param_layout = QVBoxLayout()
        self.btn_layout = QHBoxLayout()

        # Buttons
        self.btn_model = QPushButton()
        self.btn_model.setText('生成IGES模型')
        self.btn_model.clicked.connect(self.gen_model)

        self.btn_mesh = QPushButton()
        self.btn_mesh.setText('生成网格')
        self.btn_mesh.clicked.connect(self.gen_mesh)

        self.btn_exit = QPushButton()
        self.btn_exit.setText('退出')
        self.btn_exit.clicked.connect(self.close)

    def init_ui(self):
        self.setWindowTitle('BWB一体化设计')
        self.setWindowIcon(QIcon('../pic/plane_icon.jpg'))
        self.setLayout(self.global_layout)
        self.global_layout.addLayout(self.param_layout)
        self.global_layout.addLayout(self.btn_layout)

        self.btn_layout.addWidget(self.btn_model)
        self.btn_layout.addWidget(self.btn_mesh)
        self.btn_layout.addWidget(self.btn_exit)

        for i in range(0, self.sec_cnt):
            cur_sec = Section(self.airfoil_list[i], self.z_list[i], self.xf_list[i], self.yf_list[i], self.xt_list[i], self.yt_list[i], self.thickness_list[i])
            self.section.append(cur_sec)
            self.param_layout.addWidget(cur_sec)

    def gen_model(self):
        pass

    def gen_mesh(self):
        pass


if __name__ == '__main__':
    update_airfoil_list()
    app = QApplication(sys.argv)
    pm = BWB_UI(airfoil, z, x_front, y_front, x_tail, y_tail, thickness)
    pm.init_ui()
    pm.show()
    sys.exit(app.exec_())
