import numpy as np
from math import cos, sin, acos, sqrt
from src.nurbs.utility import equal, normalize


class Quaternion(object):
    symbol = ['', 'i', 'j', 'k']

    def __init__(self, w: float, x: float, y: float, z: float):
        """
        四元数 q = w + x*i + y*j + z*k
        单位4元数: q = w + x*i + y*j + z*k = cos(theta/2) + sin(theta/2) * u
        """

        self.comp = np.array([w, x, y, z], float)

    @property
    def w(self):
        return self.comp[0]

    @w.setter
    def w(self, _w):
        self.comp[0] = _w

    @property
    def x(self):
        return self.comp[1]

    @x.setter
    def x(self, _x):
        self.comp[1] = _x

    @property
    def y(self):
        return self.comp[2]

    @y.setter
    def y(self, _y):
        self.comp[2] = _y

    @property
    def z(self):
        return self.comp[3]

    @z.setter
    def z(self, _z):
        self.comp[3] = _z

    @property
    def real(self):
        return self.comp[0]

    @property
    def img(self):
        return self.comp[1:]

    @property
    def conj(self):
        return Quaternion(self.comp[0], -self.comp[1], -self.comp[2], -self.comp[3])

    @property
    def norm(self):
        return sqrt(sum(map(lambda t: t ** 2, self.comp)))

    @property
    def inv(self):
        """
        q^-1 = q' / |q|^2
        """

        return self.conj / sum(map(lambda t: t ** 2, self.comp))

    @classmethod
    def from_array(cls, _v):
        return cls(_v[0], _v[1], _v[2], _v[3])

    @classmethod
    def from_real_img(cls, _real, _img):
        return cls(_real, _img[0], _img[1], _img[2])

    @classmethod
    def from_3d(cls, x):
        return Quaternion(0, x[0], x[1], x[2])

    def __str__(self):
        ans = ''
        first_valid = False
        for i in range(4):
            if not equal(self.comp[i], 0.0):
                ans += ' ' if first_valid else ''
                ans += '+' if first_valid and self.comp[i] > 0 else ''
                ans += '{}{}'.format(self.comp[i], Quaternion.symbol[i])
                if not first_valid:
                    first_valid = True

        return ans

    def __eq__(self, other):
        return (self.comp == other.comp).all()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        return Quaternion.from_array(self.comp + other.comp)

    def __iadd__(self, other):
        self.comp += other.comp

    def __sub__(self, other):
        return Quaternion.from_array(self.comp - other.comp)

    def __isub__(self, other):
        self.comp -= other.comp

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            a1 = self.real
            a2 = other.real
            v1 = self.img
            v2 = other.img
            a = a1 * a2 - np.inner(v1, v2)
            t = a1 * v2 + a2 * v1 + np.cross(v1, v2)
            return Quaternion.from_real_img(a, t)
        else:
            return Quaternion.from_array(self.comp * other)

    def __imul__(self, other):
        if isinstance(other, Quaternion):
            a1 = self.real
            a2 = other.real
            v1 = self.img
            v2 = other.img
            a = a1 * a2 - np.inner(v1, v2)
            t = a1 * v2 + a2 * v1 + np.cross(v1, v2)
            self.w = a
            self.x = t[0]
            self.y = t[1]
            self.z = t[2]
        else:
            self.comp *= other

    def __truediv__(self, other):
        """
        a / b = a * b.inv
        """

        if isinstance(other, Quaternion):
            return self * other.inv
        else:
            return Quaternion.from_array(self.comp / other)

    def __itruediv__(self, other):
        if isinstance(other, Quaternion):
            self.__imul__(other.inv)
        else:
            self.comp /= other

    @property
    def isUnit(self):
        return equal(self.norm, 1.0)

    def normalize(self):
        self.comp /= self.norm

    @property
    def theta(self):
        return 2 * acos(self.comp[0])

    @property
    def u(self):
        st2 = sin(self.theta / 2)
        return np.array([self.comp[1] / st2, self.comp[2] / st2, self.comp[3] / st2])

    @classmethod
    def from_u_theta(cls, u, theta):
        t2 = theta / 2
        st = sin(t2)
        ct = cos(t2)
        nu = normalize(u)
        return cls(ct, st * nu[0], st * nu[1], st * nu[2])

    def rotate(self, x):
        """
        定义一个3维空间中的线性变换: Lq(x) = q * x * q', 其中'*'按照四元数的乘法定义，将x看成一个纯四元数，q'是q的共轭
        有3个性质：
        1: Lq(x+y) = Lq(x) + Lq(y) , Lq(a * x) = a * Lq(x) 其中x,y为3维向量，a为实数，该性质表明这是一个线性变换
        2：若q为单位4元数，则 ||Lq(x)|| = ||x||
        3：若q为单位4元数 且x平行于q的虚部, 则 Lq(x) = x 
        特别地，若q为单位4元数，Lq(x)为Rodriguez旋转公式，结果为x绕u逆时针旋转theta后的向量x'
        """

        if self.isUnit:
            ct = cos(self.theta)
            st = sin(self.theta)
            u = self.u
            return ct * x + (1 - ct) * np.dot(u, x) * u + st * np.cross(u, x)
        else:
            a = self.real
            v = self.img
            return (a ** 2 - sum(map(lambda _t: _t ** 2, v))) * x + 2 * np.dot(v, x) * v + 2 * a * np.cross(v, x)

    @property
    def rot_matrix(self):
        s, a, b, c = self.comp
        a2 = a ** 2
        b2 = b ** 2
        c2 = c ** 2
        ab = a * b
        sc = s * c
        ac = a * c
        sb = s * b
        bc = b * c
        sa = s * a

        return np.array([[1 - 2 * (b2 + c2), 2 * (ab - sc), 2 * (ac + sb)],
                         [2 * (ab + sc), 1 - 2 * (a2 + c2), 2 * (bc - sa)],
                         [2 * (ac - sb), 2 * (bc + sa), 1 - 2 * (a2 + b2)]])

    @classmethod
    def from_rot_matrix(cls, r, positive=True):
        a = 0.5 * sqrt(1 + r[0][0] + r[1][1] + r[2][2]) * (1.0 if positive else -1.0)
        b = 0.25 * (r[2][1] - r[1][2]) / a
        c = 0.25 * (r[0][2] - r[2][0]) / a
        d = 0.25 * (r[1][0] - r[0][1]) / a
        return cls(a, b, c, d)


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


class DCM(object):
    def __init__(self, base1, base2):
        """
        方向余弦矩阵
        :param base1: 起始坐标轴标架 
        :param base2: 目标坐标轴标架
        """

        I = normalize(base1[0])
        J = normalize(base1[1])
        K = normalize(base1[2])
        i = normalize(base2[0])
        j = normalize(base2[1])
        k = normalize(base2[2])

        self.dcm = np.matrix([[np.dot(I, i), np.dot(I, j), np.dot(I, k)],
                              [np.dot(J, i), np.dot(J, j), np.dot(J, k)],
                              [np.dot(K, i), np.dot(K, j), np.dot(K, k)]])

    @property
    def rot_matrix(self):
        return self.dcm
