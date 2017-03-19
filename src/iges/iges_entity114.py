from src.iges.iges_core import *

class IGES_Entity114(IGES_Entity):
    '''
    Parametric Spline Surface
    '''

    def __init__(self, _u, _v, _C):
        super(IGES_Entity114, self).__init__(114, IGES_Entity.SeqCnt + 1)

    def BuildParam(self):
        # Generate raw ASCII record without sequence number
        param=""
        param+="{},".format(self.directory.Entity_Type_Number)

        return self.ConvertRawToFormatted(param)

    def BuildSection(self):
        self.directory_record = self.directory.BuildEntry()
        self.param_record = self.BuildParam()

class IGES_Entity114_Builder(IGES_Entity_Builder):

    def __init__(self):
        pass

    def GetEntity(self):
        return self.entity