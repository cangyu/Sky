import numpy as np


class XF_Component(object):
    def __init__(self, index):
        self.index = index


class XF_Comment(XF_Component):
    def __init__(self, msg=''):
        super(XF_Comment, self).__init__(0)
        self.msg = msg


class XF_Header(XF_Component):
    def __init__(self, msg='Grid generated with Python 3.x'):
        super(XF_Header, self).__init__(1)
        self.msg = msg


class XF_Dimension(XF_Component):
    def __init__(self, dim):
        super(XF_Dimension, self).__init__(2)
        self.dim = dim


class XF_Node(XF_Component):
    def __init__(self):
        super(XF_Node, self).__init__(10)


class XF_Cell(XF_Component):
    def __init__(self):
        super(XF_Cell, self).__init__(12)


class XF_Face(XF_Component):
    def __init__(self):
        super(XF_Face, self).__init__(13)
