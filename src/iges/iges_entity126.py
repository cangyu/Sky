import numpy as np
from iges_core import IGES_Entity

class IGES_Entity126(IGES_Entity):
    '''
    NURBS Curve Entity
    '''

    def __init__(self, p, n, planar, closed, polynomial, periodic, knots, weights, ctrl_pts, sp, ep, norm, form=0):
        super(IGES_Entity126, self).__init__(126)
        self.directory.Form_Number = form

        m = n + p + 1
        if knots.shape[0] != m + 1 or weights.shape[0] != n + 1 or ctrl_pts.shape != (n + 1, 3) or norm.shape[0] != 3:
            raise Exception("Invalid description!")

        self.K = int(n)
        self.M = int(p)
        self.PROP1 = int(planar)
        self.PROP2 = int(closed)
        self.PROP3 = int(polynomial)
        self.PROP4 = int(periodic)

        self.T = np.zeros(m + 1, float)
        self.W = np.zeros(n + 1, float)
        self.X = np.zeros(n + 1, float)
        self.Y = np.zeros(n + 1, float)
        self.Z = np.zeros(n + 1, float)

        self.V = np.array([float(sp), float(ep)])
        self.XNORM = float(norm[0])
        self.YNORM = float(norm[1])
        self.ZNORM = float(norm[2])

        for i in range(0, m + 1):
            self.T[i] = knots[i]

        for i in range(0, n + 1):
            self.W[i] = weights[i]
            self.X[i] = ctrl_pts[i][0]
            self.Y[i] = ctrl_pts[i][1]
            self.Z[i] = ctrl_pts[i][2]

    def BuildParam(self):
        param = ""
        param += "{},".format(self.directory.Entity_Type_Number)
        param += "{},{},".format(self.K, self.M)
        param += "{},{},{},{},".format(self.PROP1, self.PROP2, self.PROP3, self.PROP4)

        for i in range(0, len(self.T)):
            param += "{},".format(self.T[i])

        for i in range(0, len(self.W)):
            param += "{},".format(self.W[i])

        for i in range(0, len(self.X)):
            param += "{},{},{},".format(self.X[i], self.Y[i], self.Z[i])

        param += "{},{},{},{},{};".format(self.V[0], self.V[1], self.XNORM, self.YNORM, self.ZNORM)

        return self.ConvertRawToFormatted(param)
