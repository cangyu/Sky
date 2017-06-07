import numpy as np
from abc import abstractmethod


class XF_Section(object):
    def __init__(self, index):
        self.index = index
        self.formatted_content = None

    @abstractmethod
    def build_content(self):
        pass

    @property
    def content(self):
        return self.formatted_content


class XF_Comment(XF_Section):
    def __init__(self, msg=''):
        super(XF_Comment, self).__init__(0)
        self.msg = msg

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.msg)


class XF_Header(XF_Section):
    def __init__(self, info='Grid generated with Python 3.x'):
        super(XF_Header, self).__init__(1)
        self.info = info

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.info)


class XF_Dimension(XF_Section):
    def __init__(self, ND):
        super(XF_Dimension, self).__init__(2)
        self.ND = ND

    def build_content(self):
        self.formatted_content = "({} {})".format(self.index, self.ND)


class XF_Node(XF_Section):
    def __init__(self, zone, first, last, node_type, ND, pts=None):
        super(XF_Node, self).__init__(10)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.node_type = node_type
        self.ND = ND
        self.pts = None if pts is None else np.copy(pts)

    def build_content(self):
        self.formatted_content = "({} ({} {} {} {} {})".format(self.index,
                                                               hex(self.zone_id)[2:],
                                                               hex(self.first_index)[2:],
                                                               hex(self.last_index)[2:],
                                                               hex(self.node_type)[2:],
                                                               self.ND)
        if self.pts is not None:
            n, dim = self.pts.shape
            self.formatted_content += "(\n"
            for i in range(n):
                for d in range(dim):
                    self.formatted_content += "{} ".format(self.pts[i][d])
                self.formatted_content += '\n'
            self.formatted_content += ")"
        self.formatted_content += ")"


class XF_Cell(XF_Section):
    def __init__(self, zone, first, last, cell_type, elem_type, cell_info=None):
        super(XF_Cell, self).__init__(12)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.cell_type = cell_type
        self.elem_type = elem_type
        self.cell_info = None if cell_info is None else np.copy(cell_info)

    def build_content(self):
        self.formatted_content = "({} ({} {} {} {} {})".format(self.index,
                                                               hex(self.zone_id)[2:],
                                                               hex(self.first_index)[2:],
                                                               hex(self.last_index)[2:],
                                                               self.cell_type,
                                                               self.elem_type)
        if self.cell_info is not None:
            self.formatted_content += '(\n'
            for ci in self.cell_info:
                self.formatted_content += "{} ".format(ci)
            self.formatted_content += '\n)'

        self.formatted_content += ')'


class XF_Face(XF_Section):
    def __init__(self):
        super(XF_Face, self).__init__(13)

    def build_content(self):
        pass


class XF_MSH(object):
    def __init__(self):
        self.section_list = []

    def add_section(self, section: XF_Section):
        self.section_list.append(section)

    def write(self, fn):
        msh = open(fn, 'w')
        for i in range(len(self.section_list)):
            self.section_list[i].build_content()
            msh.write("{}\n".format(self.section_list[i].content))

        msh.close()
