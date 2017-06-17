import numpy as np
from src.nurbs.utility import equal
from functools import reduce


def square(u: float):
    return u ** 2


class Quaternion(object):
    symbol = ['', 'i', 'j', 'k']

    def __init__(self, w: float, x: float, y: float, z: float):
        """
        四元数 q = w + x*i + y*j + z*k
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
        return np.sqrt(reduce(square, self.comp))

    @property
    def inv(self):
        t = self.conj
        t.comp /= reduce(square, self.comp)
        return t

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
        a1 = self.real
        a2 = other.real
        v1 = self.img
        v2 = other.img
        a = a1 * a2 - np.inner(v1, v2)
        t = a1 * v2 + a2 * v1 + np.cross(v1, v2)
        return Quaternion.from_real_img(a, t)

    def __imul__(self, other):
        a1 = self.real
        a2 = other.real
        v1 = self.img
        v2 = other.img
        a = a1 * a2 - np.inner(v1, v2)
        t = a1 * v2 + a2 * v1 + np.cross(v1, v2)
        self.comp = np.array([a, t[0], t[1], t[2]])

    def normalize(self):
        """
        单位化 
        """

        self.comp /= reduce(square, self.comp)

    def linear_transform(self, x):
        """
        定义一个3维空间中的线性变换: Lq(x) = q * x * q'
        其中'*'按照四元数的乘法定义，将x看成一个纯四元数，q'是q的共轭
        :param x: 待变换空间向量
        :return: 变换后的空间向量
        """

        a = self.real
        v = self.img
        t1 = (a ** 2 - reduce(square, v)) * x
        t2 = 2 * np.dot(v, x) * v
        t3 = 2 * a * np.cross(v, x)
        return t1 + t2 + t3
