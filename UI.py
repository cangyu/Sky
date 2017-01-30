import sys
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QLineEdit
from PyQt5.QtWidgets import QGridLayout
from PyQt5.QtWidgets import QPushButton

class ParametricModeling(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self)
        self.setWindowTitle('简单客机一体化设计')

        tjtw=QLabel('机头宽度(mm):')
        ejtw=QLineEdit()

        tjth=QLabel('机头宽度(mm):')
        ejth = QLineEdit()

        tjtl=QLabel('机头长度(mm):')
        ejtl = QLineEdit()

        btn_gen=QPushButton(self)
        btn_gen.setText('生成')
        btn_gen.setToolTip('根据输入参数生成模型')


        btn_exit=QPushButton(self)
        btn_exit.setText('退出')
        btn_exit.setToolTip('退出本程序')

        grid=QGridLayout()
        grid.setSpacing(10)

        grid.addWidget(tjtl, 0, 0)
        grid.addWidget(ejtl, 0, 1)
        grid.addWidget(tjtw, 1, 0)
        grid.addWidget(ejtw, 1, 1)
        grid.addWidget(tjth, 2, 0)
        grid.addWidget(ejth, 2, 1)
        grid.addWidget(btn_gen, 3, 0)
        grid.addWidget(btn_exit, 3, 1)

        self.setLayout(grid)
        self.resize(350, 300)


if __name__ == '__main__':
    app=QApplication(sys.argv)

    ythsj=ParametricModeling()
    ythsj.show()

    sys.exit(app.exec_())

