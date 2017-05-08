from iges_core import *

class IGES_Entity116(IGES_Entity):
    '''
    Point Entity
    '''

    def __init__(self, _x, _y, _z, _ptr=0):
        super(IGES_Entity116, self).__init__(116)

        self.X = float(_x)
        self.Y = float(_y)
        self.Z = float(_z)
        self.PTR = int(_ptr)

    def BuildParam(self):
        param=""
        param+="{},".format(self.directory.Entity_Type_Number)
        param+="{},{},{},{};".format(self.X,self.Y,self.Z,self.PTR)

        return self.ConvertRawToFormatted(param)
