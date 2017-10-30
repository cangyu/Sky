import numpy as np

from nurbs import point_inverse
from src.aircraft.frame import HWBFrame

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

'''Wing Planform'''
cr = 19.2
spn = 21
fl = np.array([2.4, 3.0, 3.0])
alpha = np.radians([35, 60, 36, 28])
tl = np.array([1.5, 4.5, 2.5])
beta = np.radians([-27, -42, -10, 12.5])
frm = HWBFrame(cr, spn, fl, alpha, tl, beta)

fc = frm.front_crv
tc = frm.tail_crv

p = fc(0.23)

u = point_inverse(tc, p[2], 2)

print(p)
print(u)
print(tc(u))
