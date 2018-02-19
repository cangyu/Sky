import numpy as np
from matplotlib import pyplot as plt
from planform import HWBNoseBluntPlanform
from spacing import chebshev_dist_multi


spn2 = 21.0
root_chord = 19.2
middle_chord = 4.50
tip_chord = 1.9
leading_cpt = [(0.18, 0.8), (1.97, 2.2), (5.12, 4.17), (7.14, 5.47),
               (8.58, 6.95), (10.2, 9.6), (17.4, spn2)]
trailing_cpt = [(15.15, 8.25), (15.1, 10.85)]

_chord = [root_chord, middle_chord, tip_chord]
_cpt = leading_cpt + trailing_cpt
planform = HWBNoseBluntPlanform(_chord, _cpt)

u = [0.0]
for p in _cpt:
    r = p[1]/spn2
    if r not in u:
        u.append(r)

u.sort()

_seg = [0.] + [leading_cpt[i][1] / spn2 for i in [0, 1, 2, 5, 6]]
_num = [3, 3, 5, 5, 7]
# u = chebshev_dist_multi(_seg, _num)
print(u)
# u = np.array([0, 2.143, 4.285, 6.428, 8.571, 0.09, 0.18, 0.26, 0.38, 0.51, 0.62, 0.74, 0.85, 0.93, 1])
n = len(u)
z = np.array([_r * spn2 for _r in u])

fig = plt.figure()
ax = fig.add_subplot(111)
planform.pic(ax, u=u)
fig.tight_layout()
plt.show()
