from src.iges.iges_core import *


class IGES_Entity110(IGES_Entity):
    '''
    Line Entity
    '''

    def __init__(self, _p1, _p2, _form=0):
        super(IGES_Entity110, self).__init__(110, IGES_Entity.SeqCnt + 1)
        self.directory.Form_Number = _form
        self.X1 = float(_p1[0])
        self.Y1 = float(_p1[1])
        self.Z1 = float(_p1[2])
        self.X2 = float(_p2[0])
        self.Y2 = float(_p2[1])
        self.Z2 = float(_p2[2])

    def BuildParam(self):
        # Generate raw ASCII record without sequence number
        param = ""
        param += "{},".format(self.directory.Entity_Type_Number)
        param += "{},{},{},".format(self.X1, self.Y1, self.Z1)
        param += "{},{},{};".format(self.X2, self.Y2, self.Z2)

        return self.ConvertRawToFormatted(param)

    def BuildSection(self):
        self.directory_record = self.directory.BuildEntry()
        self.param_record = self.BuildParam()


class IGES_Entity110_Builder(IGES_Entity_Builder):
    '''
    Line Builder
    '''

    def __init__(self, _ends, _f=0):
        self.entity = IGES_Entity110(_ends[0], _ends[1], _f)
        self.entity.BuildSection()

    def GetEntity(self):
        return self.entity
