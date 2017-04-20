import os
import sys
import win32com.client

"""
documents1 = app.Documents
partDocument1 = documents1.Add("Part")
part1 = partDocument1.Part
bodies1 = part1.Bodies
body1 = bodies1.Item("PartBody")
sketches1 = body1.Sketches
originElements1 = part1.OriginElements
reference1 = originElements1.PlaneXY
sketch1 = sketches1.Add(reference1)
"""


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
