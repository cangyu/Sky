from PyQt5.QtWidgets import *
from aeroplane import Aircraft
import sys


class ParametricModeling(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self)

        #aircraft designing
        self.aircraft=Aircraft()

        # buttons
        self.btn_calc = QPushButton(self)
        self.btn_calc.setText('计算')
        self.btn_calc.setToolTip('根据原始参数计算衍生参数')
        self.btn_calc.clicked.connect(self.update_param)

        self.btn_gen = QPushButton(self)
        self.btn_gen.setText('生成')
        self.btn_gen.setToolTip('根据输入参数生成模型')
        self.btn_gen.clicked.connect(self.gen_model)

        self.btn_exit = QPushButton(self)
        self.btn_exit.setText('退出')
        self.btn_exit.setToolTip('退出本程序')
        self.btn_exit.clicked.connect(self.close)

        # layout
        self.layout =QVBoxLayout()
        self.param_layout = QVBoxLayout()
        self.btn_layout = QHBoxLayout()

        self.init_UI()

    def init_UI(self):

        self.aircraft.init_widget()

        self.setWindowTitle('简单客机一体化设计')
        self.setLayout(self.layout)

        self.layout.addLayout(self.param_layout)
        self.layout.addLayout(self.btn_layout)

        self.param_layout.addWidget(self.aircraft.widget)

        self.btn_layout.addWidget(self.btn_calc)
        self.btn_layout.addWidget(self.btn_gen)
        self.btn_layout.addWidget(self.btn_exit)

    def update_param(self):
        self.aircraft.update_derived_param()

    def gen_model(self):
        self.aircraft.generate()

if __name__ == '__main__':
    app = QApplication(sys.argv)

    ythsj = ParametricModeling()
    ythsj.show()

    sys.exit(app.exec_())
