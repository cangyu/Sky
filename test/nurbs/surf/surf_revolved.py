import numpy as np
from src.nurbs.curve import Arc, Line
from src.nurbs.surface import RevolvedSurf
from src.iges.iges_core import IGES_Model

try:
    from src.misc.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True

fn = "revolved.igs"
model_file = IGES_Model(fn)

generatrix1 = Line((50, 0, 0), (50, 0, 20))
model_file.add_entity(generatrix1.to_iges())

generatrix2 = Arc.from_2pnt((100, 0, 0), (100, 0, 20), 180, (1, 0, 0))
model_file.add_entity(generatrix2.to_iges())

revolved1 = RevolvedSurf((20, 0, 0), (0, 0, 1), 60, generatrix1)
model_file.add_entity(revolved1.to_iges())

revolved2 = RevolvedSurf((50, 0, 0), (0, 0, 1), 180, generatrix2)
model_file.add_entity(revolved2.to_iges())

model_file.write()

if auto_view:
    view(fn)
