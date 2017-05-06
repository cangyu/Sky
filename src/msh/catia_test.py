import win32com.client.dynamic
import sys
import os
import numpy as np
import win32gui

halfSpan = 1000.0
footLength = 100.0
tipLength = 50.0
rootTwist = 0.0
tipTwist = 5.0

RootAirfoil = np.loadtxt('../../airfoil/NACA0012.dat')
TipAirfoil = np.loadtxt('../../airfoil/NACA0012.dat')

CATIA = win32com.client.Dispatch("CATIA.Application")
BWB_DOC = CATIA.Documents
PART_DOC = BWB_DOC.Add("Part")

part = PART_DOC.Part
ShFactory = part.HybridShapeFactory

bodies1 = part.HybridBodies
body = bodies1.Add()
body.name = "WireFrame"

bodies2 = body.hybridBodies
body2 = bodies2.Add()
body2.Name = "RootSection"
body3 = bodies2.Add()
body3.Name = "TipSection"
body4 = bodies1.Add()
body4.Name = "Surface"

point0 = ShFactory.AddNewPointCoord(0.0, 0.0, 0.0)
body.AppendHybridShape(point0)
part.Update()

wing_axis = ShFactory.AddNewDirectionByCoord(0.0, 0.0, 1.0)
