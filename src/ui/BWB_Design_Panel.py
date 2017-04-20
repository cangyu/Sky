import sys
import os
import math
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from src.aircraft.wing import *
from src.aircraft.aircraft import *
from src.com.catia import *


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
        self.airfoil_combobox.addItems(AIRFOIL_LIST)
        self.airfoil_combobox.setCurrentText(self.airfoil)
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

    def update(self):
        self.airfoil = self.airfoil_combobox.currentText()
        self.thickness = float(self.thick_le.text())
        self.z = float(self.z_le.text())
        self.xf = float(self.xf_le.text())
        self.yf = float(self.yf_le.text())
        self.xt = float(self.xt_le.text())
        self.yt = float(self.yt_le.text())


class BWB_UI(QDialog):
    def __init__(self):
        super(BWB_UI, self).__init__()

        # Contents
        self.z_list = np.array([0., 0.624029, 1.38967, 2.43503, 3.73439, 5.25574, 6.96162, 8.81003, 10.7555, 12.75, 14.7445, 16.69, 18.5384, 20.2443, 21.7656, 23.065, 24.1103, 24.876, 25.343, 25.5])
        self.xf_list = np.array([0., 0.05, 0.3, 1.7, 4.9, 6.85, 8.45, 9.65, 10.6, 11.1, 11.7, 12.1, 12.4, 12.8, 13.2, 13.7, 14.1, 14.5, 15.2, 16.])
        self.yf_list = np.zeros(len(self.z_list), float)
        self.xt_list = np.array([19.7, 19.6, 19.6, 19.5, 19.3, 19, 18.3, 17.3, 16.6, 16.5, 16.8, 17, 17.45, 17.8, 18.1, 18.4, 18.55, 18.65, 18.3, 17.8])
        self.yt_list = np.zeros(len(self.z_list), float)
        self.thickness_list = np.ones(len(self.z_list), float)
        self.row_cnt = len(self.z_list)
        self.col_cnt = len(BWB_Section.BWB_SEC_PARAM)
        self.airfoil_list = []
        self.section = []

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

        # Init Content
        for i in range(0, self.row_cnt):
            csec = BWB_Section(AIRFOIL_LIST[0], self.thickness_list[i], self.z_list[i], self.xf_list[i], self.yf_list[i], self.xt_list[i], self.yt_list[i])
            self.section.append(csec)
            for j in range(0, self.col_cnt):
                self.section_param_list.setCellWidget(i, j, csec.widget_list[j])

    def get_content(self):
        self.airfoil_list = []
        self.thickness_list = []
        self.z_list = []
        self.xf_list = []
        self.yf_list = []
        self.xt_list = []
        self.yt_list = []
        for i in range(0, self.row_cnt):
            self.section[i].update()
            self.airfoil_list.append(self.section[i].airfoil)
            self.thickness_list.append(self.section[i].thickness)
            self.z_list.append(self.section[i].z)
            self.xf_list.append(self.section[i].xf)
            self.yf_list.append(self.section[i].yf)
            self.xt_list.append(self.section[i].xt)
            self.yt_list.append(self.section[i].yt)

    def gen_model(self):
        self.get_content()
        model = Aircraft(self.airfoil_list, self.thickness_list, self.z_list, self.xf_list, self.yf_list, self.xt_list, self.yt_list)
        fileout = model.generate(True, True, False, False)
        view(fileout)

    def gen_mesh(self):
        pass


if __name__ == '__main__':
    update_airfoil_list()
    app = QApplication(sys.argv)
    pm = BWB_UI()
    pm.show()
    sys.exit(app.exec_())
