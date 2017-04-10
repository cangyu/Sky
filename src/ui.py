from PyQt5.QtWidgets import *
from PyQt5.QtGui import QIcon
import sys
import os
import math

# 翼型文件夹路径
airfoil_dir = '../airfoil/'

# 翼型列表
airfoil_list = []


def update_airfoil_list():
    for f in os.listdir(airfoil_dir):
        cur_filename = os.path.join(airfoil_dir, f)
        if os.path.isfile(cur_filename) and cur_filename.endswith('.dat'):
            airfoil_list.append(f)


class AircraftUI(QWidget):
    def __init__(self):
        QWidget.__init__(self)

    def initUI(self):
        self.setWindowTitle('BWB一体化设计')
        self.setWindowIcon(QIcon('../pic/plane_icon.jpg'))


if __name__ == '__main__':
    app = QApplication(sys.argv)

    update_airfoil_list()

    pm = AircraftUI()
    pm.initUI()
    pm.show()

    sys.exit(app.exec_())
