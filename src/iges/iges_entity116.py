from src.iges.iges_core import  IGES_Entity

class IGES_Entity116(IGES_Entity):
    '''
    Point Entity
    '''

    def __init__(self, _x, _y, _z, _ptr=0):
        super().__init__(116)

        self.X = float(_x)
        self.Y = float(_y)
        self.Z = float(_z)
        self.PTR = int(_ptr)

    def toAsciiParam(self):
        pass