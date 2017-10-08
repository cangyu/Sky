from src.iges.iges_core import Model
from src.iges.iges_entity110 import Entity110

a = [0, 0, 0]
b = [10, 20, 30]

line = Entity110(a, b)
model = Model()
model.add_entity(line)
model.save('line.igs')
