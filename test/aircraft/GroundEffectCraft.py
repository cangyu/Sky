import numpy as np
import math
from src.iges.iges_core import IGES_Model

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

fn = "GroundEffectCraft.igs"
gec = IGES_Model(fn)

fuselage_width = 2.93
fuselage_height = 3.09

hu = 2.472
hd = hu - fuselage_height
nl = 4.66
nfr = 0.05
r1 = 0.9
theta1 = math.radians(15)
l1 = 0.66

p0 = np.array([0, nfr, 0])
p1 = np.array([r1 * (1 - math.cos(theta1)), nfr + r1 * math.sin(theta1), 0])
p2 = np.array([nl - l1, hu, 0])
p3 = np.array([nl, hu, 0])

r3 = 0.26
theta2 = math.radians(27)
l2 = 1.66

p4 = np.array([0, -nfr, 0])
p5 = np.array([r3 * (1 - math.cos(theta2)), -(nfr + r3 * math.sin(theta2)), 0])
p6 = np.array([nl - l2, hd, 0])
p7 = np.array([nl, hd, 0])

p8 = np.array([0, 0, nfr])
t8 = np.array([0, 0, 1.0])
p9 = 0.5 * (p3 + p7)
p9[-1] = fuselage_width / 2


