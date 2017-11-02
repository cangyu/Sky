from nurbs import Arc, Line
from nurbs import Coons
from src.aircraft.wing import WingProfile
from src.iges import IGES_Model

try:
    from src.misc.catia import view
except ImportError:
    auto_view = False
    print('Win32 required for CATIA usage!')
else:
    auto_view = True

u0 = WingProfile('M6', [[0, 0, 5], [20, 0, 5]], 1.6, 5).nurbs_rep
u1 = Arc.from_2pnt([25, 60, 5], [25, -60, 5], 180, [0, 0, 1])
v0 = Line(u0.start, u1.start)
v1 = Line(u0.end, u1.end)
s = Coons(u0, u1, v0, v1)

print(s)

fn = 'coons_test.igs'
model_file = IGES_Model(fn)
model_file.add_entity(s.to_iges())
model_file.add_entity(u0.to_iges())
model_file.add_entity(u1.to_iges())
model_file.add_entity(v0.to_iges())
model_file.add_entity(v1.to_iges())
model_file.write()

if auto_view:
    view(fn)
