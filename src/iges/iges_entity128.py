from src.iges.iges_core import IGES_Entity
import numpy as np


class IGES_Entity128(IGES_Entity):
    '''
    NURBS Surface Entity
    '''

    def __init__(self):
        super(IGES_Entity128, self).__init__(128)
