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

        # Add sequence number and pointer back to directory
        fp = ""
        tl = len(param)
        cs = 0
        cc = 0

        while (tl):
            cc += 1
            cl = min(64, tl)
            ce = cs + cl
            fp += "{:64} {:7}P{:7}\n".format(param[cs:ce], self.directory.Sequence_Number, IGES_Entity.SeqCnt + cc)
            tl -= cl
            cs += cl

        # Update Entity param section sequence counter
        IGES_Entity.SeqCnt += cc

        return fp

    def BuildSection(self):
        self.directory_record = self.directory.BuildEntry()
        self.param_record = self.BuildParam()

class IGES_Entity116_Builder(IGES_Entity_Builder):

    def __init__(self, _x, _y, _z):
        self.entity=IGES_Entity116(_x, _y, _z)
        self.entity.BuildSection()

    def GetEntity(self):
        return self.entity
