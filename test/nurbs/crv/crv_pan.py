from src.iges import IGES_Model
from src.geom.curve import Arc

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

dir1 = [0, 0, 300]
dir2 = [200, 0, 0]

arc1 = Arc(200, 180)
arc2 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
arc3 = Arc(200, 180)
arc3.pan(dir1)
arc4 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
arc4.pan(dir2)

model = IGES_Model('test.igs')
model.add_entity(arc1.to_iges(1, 0, [0, 0, 1]))
model.add_entity(arc2.to_iges(1, 0, [0, 1, 0]))
model.add_entity(arc3.to_iges(1, 0, [0, 0, 1]))
model.add_entity(arc4.to_iges(1, 0, [0, 1, 0]))
model.write()

if auto_view:
    view('test.igs')
