import os
import sys
import win32com.client


def view(file):
    try:
        app = win32com.client.Dispatch('catia.application')
    except:
        print("Failed to open CATIA.")
        sys.exit(-1)

    app.Visible = True
    app.DisplayFileAlerts = True

    doc = app.Documents.Open(os.path.abspath(file))


def convert(file, target_format):
    dir, filename = os.path.split(file)
    base, ext = os.path.splitext(filename)
    fileout = os.path.join(dir, base + '.' + target_format)

    try:
        app = win32com.client.Dispatch('catia.application')
    except:
        print("Failed to open CATIA.")
        sys.exit(-1)

    app.DisplayFileAlerts = False
    doc = app.Documents.Open(os.path.abspath(file))
    try:
        doc.ExportData(fileout, format)
    except:
        print("Failed to convert file: {}.".format(file))
        sys.exit(-2)
    finally:
        doc.Close()
        app.Quit()
        print("Conversion Done!")

    return fileout
