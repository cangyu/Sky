from src.nurbs.curve import Arc
from src.iges.iges_core import IGES_Model

model = IGES_Model('test.igs')

arc = Arc.from_2pnt([0, 0, 1000], [10, 0, 990], 90, [0, 1, 0])

model.add_entity(arc.to_iges(1, 0, [0, 0, 1]))

model.write()
