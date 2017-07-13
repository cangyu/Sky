from src.aircraft.wing import WingProfile
from src.nurbs.surface import RuledSurf
from src.iges.iges_core import IGES_Model

try:
    from src.misc.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True

foil1 = WingProfile('M6', [[0, 0, 0], [10, 0, 0]], p=4).nurbs_rep
foil2 = WingProfile('NACA0012', [[0, 0, 80], [15, 0, 80]], p=3).nurbs_rep
surf = RuledSurf(foil1, foil2)

fn = 'ruled.igs'
model_file = IGES_Model(fn)
model_file.add_entity(surf.to_iges(0, 0, 0, 0))
model_file.write()

if auto_view:
    view(fn)
