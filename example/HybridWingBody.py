import numpy as np
from matplotlib import pyplot as plt
from planform import HWBCommonPlanform
from airfoil import Airfoil
from profile import ProfileList
from iges import IGES_Model
from load_dist import calc_global_cl

spn2 = 21.0
root_chord = 19.2
middle_chord = 4.50
tip_chord = 2.0
leading_cpt = [(0.18, 0.8), (1.97, 2.2),
               (5.12, 4.17), (7.14, 5.47),
               (8.58, 6.95), (10.2, 9.6), (17.4, spn2)]
trailing_cpt = [(17.23, 3.83), (16.1, 6.97), (16.3, 10.15)]

_chord = [root_chord, tip_chord]
_cpt = leading_cpt + trailing_cpt
planform = HWBCommonPlanform(_chord, _cpt)

payload = 160
rho = 0.4135  # H=10000m
a = 299.5
ma = 0.8
v = ma * a
re = 4e7
area = planform.area * 0.8
cl = calc_global_cl(payload, area, ma, rho, a)
print('Design Lift Coefficient: {:.2f}'.format(cl))

pos = [0, 1.91, 3.81, 5.92, 8.03,
       10.48, 14.52, 18.24, 22.86, 26.05,
       29.85, 33.10, 36.91, 40.72, 44.56,
       48.33, 54.79, 61.25, 67.71, 74.17,
       80.62, 87.08, 93.54, 97.25, 100]
u = np.array([0.01 * x for x in pos])
n = len(u)
z = np.array([x * spn2 for x in u])

fig = plt.figure()
ax = fig.add_subplot(111)
planform.pic(ax, u=u)
fig.tight_layout()
plt.show()

sc0712 = Airfoil.from_local('SC(2)-0712')
sc0412 = Airfoil.from_local('SC(2)-0412')
sc0012 = Airfoil.from_local('SC(2)-0012')

foil = [sc0712] * n
twist = np.zeros(n)
twist_ref = np.ones(n)
wpl = ProfileList.from_planform(planform, foil, twist, twist_ref, u)

model = IGES_Model()
pfl = wpl.generate_nurbs_crv_list()
for c in pfl:
    model.add(c.to_iges())
model.save('HWB.igs')
