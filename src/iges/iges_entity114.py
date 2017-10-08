import numpy as np
from ..iges.iges_core import Entity


class Entity114(Entity):
    def __init__(self, _u, _v, _c):
        """
        Parametric Spline Surface
        :param _u: Knot vector in U direction.
        :param _v: Knot vector in V direction.
        :param _c: Coefficients.
        """

        super(Entity114, self).__init__(114)

        self.CTYPE = int(3)
        self.PTYPE = 0
        self.M = len(_u) - 1
        self.N = len(_v) - 1

        self.Tu = np.zeros(len(_u))
        for i in range(0, len(_u)):
            self.Tu[i] = _u[i]

        self.Tv = np.zeros(len(_v))
        for j in range(0, len(_v)):
            self.Tv[i] = _v[i]

        self.Coef = np.zeros((self.M + 1, self.N + 1, 3, 4, 4))

        assert _c.shape == self.Coef.shape

        for m in range(0, self.M + 1):
            for n in range(0, self.N + 1):
                for dim in range(0, 3):
                    for i in range(0, 4):
                        for j in range(0, 4):
                            self.Coef[m][n][dim][i][j] = _c[m][n][dim][i][j]

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        param = "{},".format(self.directory.entity_type_number)
        param += "{},".format(self.CTYPE)
        param += "{},".format(self.PTYPE)
        param += "{},".format(self.M)
        param += "{},".format(self.N)

        for i in range(0, self.M + 1):
            param += "{},".format(self.Tu[i])

        for i in range(0, self.N + 1):
            param += "{},".format(self.Tv[i])

        for m in range(0, self.M + 1):
            for n in range(0, self.N + 1):
                for dim in range(0, 3):
                    for i in range(0, 4):
                        for j in range(0, 4):
                            if m == self.M and n == self.N and dim == 2 and i == 3 and j == 3:
                                param += "{};".format(self.Coef[m][n][dim][i][j])
                            else:
                                param += "{},".format(self.Coef[m][n][dim][i][j])

        return self.to_formatted(param)
