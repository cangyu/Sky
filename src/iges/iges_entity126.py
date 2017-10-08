import numpy as np
from ..iges.iges_core import Entity


class Entity126(Entity):
    def __init__(self, p, n, planar, closed, polynomial, periodic, knots, weights, ctrl_pts, sp, ep, norm, form=0):
        """
        NURBS Curve Entity
        :param p: Degree of basis functions.
        :type p: int
        :param n: The last index of control points.
        :type n: int
        :param planar: 0 = non-planar, 1 = planar
        :type planar: int
        :param closed: 0 = open curve, 1 = closed curve
        :type closed: int
        :param polynomial: 0 = rational, 1 = polynomial
        :type polynomial: int
        :param periodic: 0 = non-periodic, 1 = periodic
        :type periodic: int
        :param knots: Knot vector.
        :param weights: Rational weight coefficients.
        :param ctrl_pts: Control points.
        :param sp: Starting point.
        :param ep: Ending point.
        :param norm: Unit normal vector. (If curve is planar)
        :param form: Form number. (0-5).
        :type form: int
        """

        super(Entity126, self).__init__(126)
        self.directory.Form_Number = form

        m = n + p + 1
        if len(knots) != m + 1:
            raise ValueError("Invalid Knot Vector!")
        if len(weights) != n + 1:
            raise ValueError("Invalid Weights!")
        if ctrl_pts.shape != (n + 1, 3):
            raise ValueError("Invalid Control Points!")
        if len(norm) != 3:
            raise ValueError("Invalid Norm!")

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

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        param = "{},".format(self.directory.entity_type_number)
        param += "{},{},".format(self.K, self.M)
        param += "{},{},{},{},".format(self.PROP1, self.PROP2, self.PROP3, self.PROP4)

        for i in range(0, len(self.T)):
            param += "{},".format(self.T[i])

        for i in range(0, len(self.W)):
            param += "{},".format(self.W[i])

        for i in range(0, len(self.X)):
            param += "{},{},{},".format(self.X[i], self.Y[i], self.Z[i])

        param += "{},{},{},{},{};".format(self.V[0], self.V[1], self.XNORM, self.YNORM, self.ZNORM)

        return self.to_formatted(param)
