import numpy as np
from math import sin, cos


class EulerAngle(object):
    def __init__(self, a, b, g):
        """
        Intrinsic Rotation
        :param a: rotation angle around z-axis
        :param b: rotation angle around y-axis
        :param g: rotation angle around x-axis
        """

        self.alpha = a
        self.beta = b
        self.gamma = g

    @property
    def roll(self):
        return self.gamma

    @roll.setter
    def roll(self, val):
        self.gamma = val

    @property
    def pitch(self):
        return self.beta

    @pitch.setter
    def pitch(self, val):
        self.beta = val

    @property
    def yaw(self):
        return self.alpha

    @yaw.setter
    def yaw(self, val):
        self.alpha = val

    @property
    def z_rot_matrix(self):
        sa = sin(self.alpha)
        ca = cos(self.alpha)

        return np.matrix([ca, -sa, 0],
                         [sa, ca, 0],
                         [0, 0, 1])

    @property
    def y_rot_matrix(self):
        sb = sin(self.beta)
        cb = cos(self.beta)

        return np.matrix([[cb, 0, sb],
                          [0, 1, 0],
                          [-sb, 0, cb]])

    @property
    def x_rot_matrix(self):
        sg = sin(self.gamma)
        cg = cos(self.gamma)

        return np.matrix([[1, 0, 0],
                          [0, cg, -sg],
                          [0, sg, cg]])

    @property
    def rot_matrix(self):
        """
        R(alpha, beta, gamma) = Rz(alpha) * Ry(beta) * Rx(gamma)
        :return: Rotation matrix
        """

        sa = sin(self.alpha)
        ca = cos(self.alpha)
        sb = sin(self.beta)
        cb = cos(self.beta)
        sg = sin(self.gamma)
        cg = cos(self.gamma)

        return np.matrix([[ca * cb, ca * sb * sg - sa * cg, ca * sb * cg + sa * sg],
                          [sa * cb, sa * sb * sg + ca * cg, sa * sb * cg - ca * sg],
                          [-sb, cb * sg, cb * cg]])
