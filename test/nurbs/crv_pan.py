from src.nurbs.curve import Arc
from src.iges.iges_core import IGES_Model

model = IGES_Model('test.igs')

dir1 = [0, 0, 300]
dir2 = [200, 0, 0]

arc1 = Arc(200, 180)
arc2 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
arc3 = Arc(200, 180)
arc3.pan(dir1)
arc4 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
arc4.pan(dir2)

model.add_entity(arc1.to_iges(1, 0, [0, 0, 1]))
model.add_entity(arc2.to_iges(1, 0, [0, 1, 0]))
model.add_entity(arc3.to_iges(1, 0, [0, 0, 1]))
model.add_entity(arc4.to_iges(1, 0, [0, 1, 0]))

model.write()