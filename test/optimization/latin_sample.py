import numpy as np
from src.opt.latin_cube import latin_sample


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


a = (0, 1, 2, 3, 4, 5)
b = ['a', 'b', 'c', 'd', 'e', 'f']
c = ('I', 'II', 'III', 'IV', 'V', 'VI')
d = np.array([1.1, 2.2, 3.3, 4.4, 5.5, 6.6])
e = (Person('ggsmd', 21), Person('sbtty', 22), Person('smdhx', 23), Person('shzyh', 24), Person('sgdb', 25), Person('tltbpa', 26))

sample = latin_sample([a, b, c, d, e])
for k, item in enumerate(sample):
    print("{} {}".format(k, repr(item)))
