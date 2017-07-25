import numpy as np
import math
from src.aircraft.wing import Wing, WingProfile
from src.iges.iges_core import IGES_Model
from src.nurbs.curve import Arc
from src.nurbs.surface import ExtrudedSurf, RuledSurf

try:
    from src.misc.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True

fn = "GroundEffectCraft.igs"
gec = IGES_Model(fn)

foil = ['NACA0012', 'M6', 'M6']
z_offset = np.array([0, 4.5, 7.7])
length = np.array([3.3, 3.3, 0.9])
sweep_back = np.array([0, 0, 25], float)
twist = np.array([5, 5, 2], float)
dihedral = np.array([0, 3, 5], float)
twist_pos = np.array([0.25, 0.25, 0.25])
y_ref = np.zeros(3)
thickness_factor = np.ones(3, float)

sect_list = []
crv_list = []

for k in range(3):
    wp = WingProfile.from_geom_param(foil[k], z_offset[k], length[k], sweep_back[k], twist[k], dihedral[k], twist_pos[k], y_ref[k], thickness_factor[k])
    sect_list.append(wp)
    crv = wp.nurbs_rep()
    crv_list.append(crv)

inner_surf = RuledSurf(crv_list[0], crv_list[1])
outer_surf = RuledSurf(crv_list[1], crv_list[2])

gec.add_entity(inner_surf.to_iges())
gec.add_entity(outer_surf.to_iges())

gec.write()
if auto_view:
    view(fn)
