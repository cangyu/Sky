import numpy as np
from src.iges.iges_core import Entity


class Entity112(Entity):
    def __init__(self, _t, _c):
        """
        Parametric Spline Curve
        :param _t: Number of segments
        :param _c: Coordinate polynomial
        """

        super(Entity112, self).__init__(112)

        # Spline Type
        self.CTYPE = int(3)

        # Degree of continuity with respect to arc length
        self.H = int(2)

        # Number of dimensions
        self.NDIM = int(3)

        # Number of segments
        self.N = len(_t) - 1

        # Break points of piecewise polynomial
        self.T = np.zeros(len(_t))
        for i in range(0, len(_t)):
            self.T[i] = _t[i]

        # Coordinate polynomial
        self.C = np.zeros((self.N, 3, 4))
        for i in range(0, self.N):
            for j in range(0, 3):
                for k in range(0, 4):
                    self.C[i][j][k] = _c[i][j][k]

        # Terminal info
        self.TPX0 = _c[self.N][0][0]  # X value
        self.TPX1 = _c[self.N][0][1]  # X first derivative
        self.TPX2 = _c[self.N][0][2]  # X second derivative/2!
        self.TPX3 = _c[self.N][0][3]  # X third derivative/3!

        self.TPY0 = _c[self.N][1][0]  # Y value
        self.TPY1 = _c[self.N][1][1]  # Y first derivative
        self.TPY2 = _c[self.N][1][2]  # Y second derivative/2!
        self.TPY3 = _c[self.N][1][3]  # Y third derivative/3!

        self.TPZ0 = _c[self.N][2][0]  # Z value
        self.TPZ1 = _c[self.N][2][1]  # Z first derivative
        self.TPZ2 = _c[self.N][2][2]  # Z second derivative/2!
        self.TPZ3 = _c[self.N][2][3]  # Z third derivative/3!

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        """Generate raw ASCII record without sequence number"""
        param = "{},".format(self.directory.entity_type_number)
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

        '''Convert to IGES formatted strings'''
        return self.to_formatted(param)
