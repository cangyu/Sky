from src.iges.iges_core import *

class IGES_Entity116(IGES_Entity):
    '''
    Point Entity
    '''

    def __init__(self, _x, _y, _z, _ptr=0):
        super(IGES_Entity116, self).__init__(116, IGES_Entity.SeqCnt + 1)

        self.X = float(_x)
        self.Y = float(_y)
        self.Z = float(_z)
        self.PTR = int(_ptr)

    def BuildParam(self):
        # Generate raw ASCII record without sequence number
        param=""
        param+="{},".format(self.directory.Entity_Type_Number)
        param+="{},{},{},{};".format(self.X,self.Y,self.Z,self.PTR)

        return self.ConvertRawToFormatted(param)

    def BuildSection(self):
        self.directory_record = self.directory.BuildEntry()
        self.param_record = self.BuildParam()

class IGES_Entity116_Builder(IGES_Entity_Builder):

    def __init__(self, _x, _y, _z):
        self.entity=IGES_Entity116(_x, _y, _z)
        self.entity.BuildSection()

    def GetEntity(self):
        return self.entity
