from tkinter import Tk
from time import sleep
from tkinter.messagebox import showwarning
import win32com.client as win32

warn = lambda app: showwarning(app, 'Exit?')
RANGE = range(3, 8)


def excel():
    app = 'Excel'
    xl = win32.gencache.EnsureDispatch('%s.Application' % app)
    ss = xl.Workbooks.Add()
    sh = ss.ActiveSheet
    xl.Visible = True

    sh.Cells(1, 1).Value = 'Python-to-%s Demo' % app

    for i in RANGE:
        sh.Cells(i, 1).Value = 'Line %d' % i
    sh.Cells(10, 1).Value = "That's all folks!"

    warn(app)
    ss.Close(False)
    xl.Application.Quit()


if __name__ == '__main__':
    Tk().withdraw()
    excel()
