import os
import numpy as np

'''Path'''
SRC_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.join(SRC_DIR, '..')
AIRFOIL_DIR = os.path.join(PROJECT_DIR, 'airfoil', 'database')
XFOIL_PATH = os.path.join(PROJECT_DIR, '3rd_party', 'XFOIL', 'xfoil')

'''Local Airfoils'''
AIRFOIL_LIST = []
for f in os.listdir(AIRFOIL_DIR):
    base, ext = os.path.splitext(f)
    if ext == '.dat' or ext == '.txt':
        AIRFOIL_LIST.append(base.upper())

'''Global Reference Frame'''
global_origin = np.array([0., 0., 0.])
x_axis_positive = np.array([1., 0, 0])
x_axis_negative = np.array([-1., 0, 0])
y_axis_positive = np.array([0, 1., 0])
y_axis_negative = np.array([0, -1., 0])
z_axis_positive = np.array([0, 0, 1.])
z_axis_negative = np.array([0, 0, -1.])
