import sys
import os
import math
import time
import numpy as np
import pylab as pl

from abc import ABCMeta, abstractmethod
from src.iges.iges_core import IGES_Model
from src.iges.iges_entity112 import *
from src.iges.iges_entity116 import *

z = np.array([0., 0.624029, 1.38967, 2.43503, 3.73439, 5.25574, 6.96162,
              8.81003, 10.7555, 12.75, 14.7445, 16.69, 18.5384, 20.2443, 21.7656,
              23.065, 24.1103, 24.876, 25.343, 25.5])

x_front = np.array([0., 0.05, 0.3, 1.7, 4.9, 6.85, 8.45, 9.65, 10.6, 11.1, 11.7, 12.1,
                    12.4, 12.8, 13.2, 13.7, 14.1, 14.5, 15.2, 16.])

y_front = np.linspace(0, 0, len(z))

x_tail = np.array([19.7, 19.6, 19.6, 19.5, 19.3, 19, 18.3, 17.3, 16.6, 16.5,
                   16.8, 17, 17.45, 17.8, 18.1, 18.4, 18.55, 18.65, 18.3, 17.8])

y_tail = np.linspace(0, 0, len(z))


plane = IGES_Model()

plane.AddPart(IGES_Entity112_Builder(z, x_front, y_front, z).GetEntity())
plane.AddPart(IGES_Entity112_Builder(z, x_tail, y_tail, z).GetEntity())

for i in range(0, len(z)):
    plane.AddPart(IGES_Entity116_Builder(x_front[i], y_front[i], z[i]).GetEntity())
    plane.AddPart(IGES_Entity116_Builder(x_tail[i], y_tail[i], z[i]).GetEntity())

plane.Generate()