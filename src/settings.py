import os

SRC_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.join(SRC_DIR, '..')
AIRFOIL_DIR = os.path.join(PROJECT_DIR, 'airfoil')

AIRFOIL_LIST = []
for f in os.listdir(AIRFOIL_DIR):
    base, ext = os.path.splitext(f)
    if ext == '.dat':
        AIRFOIL_LIST.append(base)
