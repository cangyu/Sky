import unittest
import numpy as np
from src.kriging import LHS


class Person(object):
    def __init__(self, _nm, _ag):
        self._name = _nm
        self._age = _ag

    @property
    def name(self):
        return self._name

    @property
    def age(self):
        return self._age

    def __repr__(self):
        return "My name is {}, I'm {} years old now!".format(self.name, self.age)


class LHSTestCase(unittest.TestCase):
    def test_sample(self):
        a = (0, 1, 2, 3, 4, 5)
        b = ['a', 'b', 'c', 'd', 'e', 'f']
        c = ('I', 'II', 'III', 'IV', 'V', 'VI')
        d = np.array([1.1, 2.2, 3.3, 4.4, 5.5, 6.6])
        e = (Person('ggsmd', 21), Person('sbtty', 22), Person('smdhx', 23), Person('shzyh', 24), Person('sgdb', 25), Person('tltbpa', 26))
        self.assertTrue(len(a) == len(b) == len(c) == len(d) == len(e))

        lhc = LHS([a, b, c, d, e])
        sp = lhc.sample()
        for i in range(len(sp)):
            print('\n {}'.format(sp[i]))


if __name__ == '__main__':
    unittest.main()
