import numpy as np
from iges_core import IGES_Entity

class IGES_Entity128(IGES_Entity):
    '''
    NURBS Surface Entity
    '''

    def __init__(self, u, v, p1, p2, n1, n2, ctrl_pts, weights, closed_u=0, closed_v=0, poly=1, periodic_u=0, periodic_v=0, us=0.0, ue=1.0, vs=0.0, ve=1.0, form=0):
        super(IGES_Entity128, self).__init__(128)
        self.directory.Form_Number = form

        if len(u) != 2 + p1 + n1 or len(v) != 2 + p2 + n2 or ctrl_pts.shape != (1 + n1, 1 + n2, 3) or weights.shape != (1 + n1, 1 + n2):
            raise Exception("Invalid Parameter!")

        self.K1 = int(n1)  # U方向控制点最后一个下标
        self.K2 = int(n2)  # V方向控制点最后一个下标
        self.M1 = int(p1)  # U方向的阶
        self.M2 = int(p2)  # V方向的阶
        self.PROP1 = int(closed_u)
        self.PROP2 = int(closed_v)
        self.PROP3 = int(poly)
        self.PROP4 = int(periodic_u)
        self.PROP5 = int(periodic_v)

        self.U = np.array([us, ue])
        self.V = np.array([vs, ve])

        self.S = np.zeros(len(u), float)
        for i in range(0, len(u)):
            self.S[i] = u[i]

        self.T = np.zeros(len(v), float)
        for i in range(0, len(v)):
            self.T[i] = v[i]

        self.W = np.zeros((self.K1 + 1, self.K2 + 1), float)
        self.X = np.zeros((self.K1 + 1, self.K2 + 1), float)
        self.Y = np.zeros((self.K1 + 1, self.K2 + 1), float)
        self.Z = np.zeros((self.K1 + 1, self.K2 + 1), float)

        for j in range(0, self.K2 + 1):
            for i in range(0, self.K1 + 1):
                self.W[i][j] = weights[i][j]
                self.X[i][j] = ctrl_pts[i][j][0]
                self.Y[i][j] = ctrl_pts[i][j][1]
                self.Z[i][j] = ctrl_pts[i][j][2]

    def BuildParam(self):
        param = ""
        param += "{},".format(self.directory.Entity_Type_Number)
        param += "{},{},{},{},".format(self.K1, self.K2, self.M1, self.M2)
        param += "{},{},{},{},{},".format(self.PROP1, self.PROP2, self.PROP3, self.PROP4, self.PROP5)

        for u in self.S:
            param += "{},".format(u)

        for v in self.T:
            param += "{},".format(v)

        for j in range(0, self.K2 + 1):
            for i in range(0, self.K1 + 1):
                param += "{},".format(self.W[i][j])

        for j in range(0, self.K2 + 1):
            for i in range(0, self.K1 + 1):
                param += "{},{},{},".format(self.X[i][j], self.Y[i][j], self.Z[i][j])

        param += "{},{},{},{};".format(self.U[0], self.U[1], self.V[0], self.V[1])

        return self.ConvertRawToFormatted(param)
