import unittest
import sys
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *


class BasicUI_Test(unittest.TestCase):
    def test_window(self):
        app = QApplication(sys.argv)
        window = QWidget()
        window.setWindowTitle("Basic Window Test")
        window.show()
        sys.exit(app.exec_())
