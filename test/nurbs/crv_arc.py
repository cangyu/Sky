from src.nurbs.curve import Arc
from src.iges.iges_core import IGES_Model

try:
    from src.com.catia import view
except ImportError:
    print('Win32 required for CATIA usage!')

auto_view = True

arc1 = Arc(200, 180)
arc2 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])

model = IGES_Model('test.igs')
model.add_entity(arc1.to_iges(1, 0, [0, 0, 1]))
model.add_entity(arc2.to_iges(1, 0, [0, 1, 0]))
model.write()

if auto_view:
    view('test.igs')
