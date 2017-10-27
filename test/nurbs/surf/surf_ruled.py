from src.aircraft.wing import WingProfile
from src.iges import IGES_Model
from src.nurbs.surface import RuledSurf

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
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
