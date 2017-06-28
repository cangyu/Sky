import math
import numpy as np
from bitarray import bitarray
from abc import abstractmethod

var_name = ['Cr', 'Cm', 'Ct', 'Bm', 'Bt', 'Alpha_m', 'Alpha_t']
var_range = np.array([[100, 120],
                      [30, 45],
                      [6, 10],
                      [14, 30],
                      [110, 160],
                      [30, 50],
                      [10, 30]], float)


class Chromosome(object):
    def __init__(self):
        self.body = None
        self.fitness = 0.0

