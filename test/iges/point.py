from src.iges.iges_core import Model
from src.iges.iges_entity116 import Entity116

pnt = Entity116(3.14, 2.718, 0.618)
model = Model()
model.add_entity(pnt)
model.save('pnt.igs')
model.save('test.igs')
