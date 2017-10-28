from src.iges import IGES_Model
from src.geom.curve import Arc

model = IGES_Model('original.igs')
arc1 = Arc(200, 180)
arc2 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
model.add_entity(arc1.to_iges(1, 0, [0, 0, 1]))
model.add_entity(arc2.to_iges(1, 0, [0, 1, 0]))
model.write()

model = IGES_Model('reversed.igs')
arc1 = Arc(200, 180)
arc1.reverse()
arc2 = Arc.from_2pnt([0, 0, 500], [0, 0, 100], 180, [0, 1, 0])
arc2.reverse()
model.add_entity(arc1.to_iges(1, 0, [0, 0, 1]))
model.add_entity(arc2.to_iges(1, 0, [0, 1, 0]))
model.write()
