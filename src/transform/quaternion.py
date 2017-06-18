import numpy as np
from math import cos, sin, acos, sqrt
from src.nurbs.utility import equal


class Quaternion(object):
    symbol = ['', 'i', 'j', 'k']

    def __init__(self, w: float, x: float, y: float, z: float):
        """
        四元数 q = w + x*i + y*j + z*k
        单位4元数: q = w + xi + yj + zk = cos(theta/2) + sin(theta/2) * u
        """

        self.comp = np.array([w, x, y, z], float)
        self.isUnit = self.is_unit()

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
        return sqrt(sum(map(lambda _t: _t ** 2, self.comp)))

    @property
    def inv(self):
        """
        q^-1 = q' / |q|^2
        """

        return self.conj / sum(map(lambda _t: _t ** 2, self.comp))

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
        nv = self.comp + other.comp
        return Quaternion.from_array(nv)

    def __iadd__(self, other):
        self.comp += other.comp

    def __sub__(self, other):
        nv = self.comp - other.comp
        return Quaternion.from_array(nv)

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
        elif isinstance(other, (float, int)):
            return Quaternion.from_array(self.comp * other)
        else:
            raise TypeError('Invalid type!')

    def __imul__(self, other):
        if isinstance(other, Quaternion):
            a1 = self.real
            a2 = other.real
            v1 = self.img
            v2 = other.img
            a = a1 * a2 - np.inner(v1, v2)
            t = a1 * v2 + a2 * v1 + np.cross(v1, v2)
            self.__init__(a, t[0], t[1], t[2])
        elif isinstance(other, (float, int)):
            tmp = self.comp * other
            self.__init__(tmp[0], tmp[1], tmp[2], tmp[3])
        else:
            raise TypeError('Invalid type!')

    def __truediv__(self, other):
        """
        a / b = a * b.inv
        """

        if isinstance(other, Quaternion):
            return self * other.inv
        elif isinstance(other, (float, int)):
            return Quaternion.from_array(self.comp / other)
        else:
            raise TypeError('Invalid type!')

    def __itruediv__(self, other):
        if isinstance(other, Quaternion):
            tmp = self * other.inv
            self.__init__(tmp.comp[0], tmp.comp[1], tmp.comp[2], tmp.comp[3])
        elif isinstance(other, (float, int)):
            tmp = self.comp / other
            self.__init__(tmp[0], tmp[1], tmp[2], tmp[3])
        else:
            raise TypeError('Invalid type!')

    def is_unit(self):
        return equal(self.norm, 1.0)

    def normalize(self):
        self.comp /= self.norm
        self.isUnit = True

    @property
    def theta(self):
        return 2 * acos(self.comp[0])

    @property
    def u(self):
        st2 = sin(self.theta / 2)
        return np.array([self.comp[1] / st2, self.comp[2] / st2, self.comp[3] / st2])

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
