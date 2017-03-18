import numpy as np
from src.iges.iges_core import  IGES_Entity

class IGES_Entity112(IGES_Entity):
    '''
    Parametric Spline Curve
    '''

    def __init__(self, T, C):
        super(IGES_Entity112, self).__init__(112, IGES_Entity.SeqCnt + 1)

        # Spline Type
        self.CTYPE = int(3)

        # Degree of continuity with respect to arc length
        self.H = int(2)

        # Number of dimensions
        self.NDIM = int(3)

        # Number of segments
        self.N = len(T) - 1

        # Break points of piecewise polynomial
        self.T = np.zeros(len(T))
        for i in range(0, len(T)):
            self.T[i] = T[i]

        # Coordinate polynomial
        self.C = np.zeros((self.N, 3, 4))
        for i in range(0, self.N):
            for j in range(0, 3):
                for k in range(0, 4):
                    self.C[i][j][k] = C[i][j][k]

        # Terminal info
        self.TPX0 = C[self.N][0][0]  # X value
        self.TPX1 = C[self.N][0][1]  # X first derivative
        self.TPX2 = C[self.N][0][2]  # X second derivative/2!
        self.TPX3 = C[self.N][0][3]  # X third derivative/3!

        self.TPY0 = C[self.N][1][0]  # Y value
        self.TPY1 = C[self.N][1][1]  # Y first derivative
        self.TPY2 = C[self.N][1][2]  # Y second derivative/2!
        self.TPY3 = C[self.N][1][3]  # Y third derivative/3!

        self.TPZ0 = C[self.N][2][0]  # Z value
        self.TPZ1 = C[self.N][2][1]  # Z first derivative
        self.TPZ2 = C[self.N][2][2]  # Z second derivative/2!
        self.TPZ3 = C[self.N][2][3]  # Z third derivative/3!

    def BuildParam(self):

        # Generate raw ASCII record without sequence number
        param = ""
        param += "{},".format(self.directory.Entity_Type_Number)
        param += "{},".format(self.CTYPE)
        param += "{},".format(self.H)
        param += "{},".format(self.NDIM)
        param += "{},".format(self.N)

        for i in range(0, len(self.T)):
            param += "{},".format(self.T[i])

        for i in range(0, self.N):
            for j in range(0, 3):
                for k in range(0, 4):
                    param += "{},".format(self.C[i][j][k])

        param += "{},".format(self.TPX0)
        param += "{},".format(self.TPX1)
        param += "{},".format(self.TPX2)
        param += "{},".format(self.TPX3)
        param += "{},".format(self.TPY0)
        param += "{},".format(self.TPY1)
        param += "{},".format(self.TPY2)
        param += "{},".format(self.TPY3)
        param += "{},".format(self.TPZ0)
        param += "{},".format(self.TPZ1)
        param += "{},".format(self.TPZ2)
        param += "{};".format(self.TPZ3)

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