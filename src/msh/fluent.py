import numpy as np
from abc import abstractmethod


class XF_Section(object):
    def __init__(self, index):
        """
        Fluent MSH文件中每个Section的抽象类
        :param index: 该Section的类别
        """

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
        """
        注释
        :param msg: 注释信息 
        """

        super(XF_Comment, self).__init__(0)
        self.msg = msg

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.msg)


class XF_Header(XF_Section):
    def __init__(self, info='Python 3.x'):
        """
        To identify the program that wrote the file.
        :param info: Header info.
        """

        super(XF_Header, self).__init__(1)
        self.info = info

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.info)


class XF_Dimension(XF_Section):
    def __init__(self, dim: int):
        """
        Dimension info of the grid.
        :param dim: Dimension 
        """

        super(XF_Dimension, self).__init__(2)
        self.ND = dim

    def build_content(self):
        self.formatted_content = "({} {})".format(self.index, self.ND)


class XF_Node(XF_Section):
    def __init__(self, zone, first, last, node_type, ND, pts=None):
        """
        网格点
        :param zone: 区域号 
        :param first: 网格点起始序号(Starting from 1)
        :param last:  网格点终止序号
        :param node_type: 网格点类型(TGrid usage)
            0 - Virtual nodes
            1 - Any type
            2 - Boundary nodes
        :param ND: Indicate the dimensionality of the node data.
        :param pts: 网格点数组，(last-first+1)个元素
        """

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
                for d in range(self.ND):
                    self.formatted_content += "{} ".format(self.pts[i][d])
                self.formatted_content += '\n'
            self.formatted_content += ")"
        self.formatted_content += ")"

    @classmethod
    def declaration(cls, num, dim):
        return cls(0, 1, num, 0, dim)


class XF_Face(XF_Section):
    def __init__(self, zone, first, last, face_type, elem_type, face_info=None):
        """
        边界
        :param zone: 区域号
        :param first: 起始序号(Starting from 1)
        :param last: 终止序号
        :param face_type: 边界属性(十进制)
            2  - interior
            3  - wall
            4  - pressure-inlet, inlet-vent, intake-fan
            5  - pressure-outlet, outlet-vent, exhaust-fan
            7  - symmetry
            8  - periodic-shadow
            9  - pressure-far-field
            10 - velocity-inlet
            12 - periodic
            14 - fan, porous-jump, radiator
            20 - mass-flow-inlet
            24 - interface
            31 - parent(hanging node)
            36 - outflow
            37 - axis
        :param elem_type: 形状类别
            0 - mixed
            2 - linear
            3 - triangular
            4 - quadrilateral
        :param face_info: 边界信息，一维数组(last-first+1)个元素，每个元素中包含了邻接关系
        """

        super(XF_Face, self).__init__(13)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.face_type = face_type
        self.elem_type = elem_type
        self.face_info = None if face_info is None else np.copy(face_info)

    def build_content(self):
        if self.face_info is None:
            self.formatted_content = "({} ({} {} {} {})".format(self.index,
                                                                hex(self.zone_id)[2:],
                                                                hex(self.first_index)[2:],
                                                                hex(self.last_index)[2:],
                                                                hex(self.elem_type)[2:])
        else:
            self.formatted_content = "({} ({} {} {} {} {})".format(self.index,
                                                                   hex(self.zone_id)[2:],
                                                                   hex(self.first_index)[2:],
                                                                   hex(self.last_index)[2:],
                                                                   hex(self.face_type)[2:],
                                                                   hex(self.elem_type)[2:])
            self.formatted_content += '(\n'
            for fc in self.face_info:
                for cfi in fc:  # current face info
                    self.formatted_content += "{} ".format(cfi)
                self.formatted_content += '\n'
            self.formatted_content += ')'

        self.formatted_content += ')'

    @classmethod
    def declaration(cls, num):
        return cls(0, 1, num, 0, 0)


class XF_Cell(XF_Section):
    def __init__(self, zone, first, last, cell_type, elem_type, cell_info=None):
        """
        单元
        :param zone: 区域号
        :param first: 起始序号(Starting from 1)
        :param last: 终止序号
        :param cell_type: Active or not.
            0  - Dead zone.
            1  - Active zone(fluid or solid).
            32 - Inactive zone.
        :param elem_type: Type of all cell.
            0 - mixed
            1 - triangular
            2 - tetrahedral
            3 - quadrilateral
            4 - hexahedral
            5 - pyramid
            6 - wedge
        :param cell_info: 单元类型信息，一维数组，(last-first+1)个元素
        """

        super(XF_Cell, self).__init__(12)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.cell_type = cell_type
        self.elem_type = elem_type
        self.cell_info = None if cell_info is None else np.copy(cell_info)

    def build_content(self):
        if self.cell_info is not None:
            self.formatted_content = "({} ({} {} {} {} {})".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], self.cell_type, self.elem_type)
            self.formatted_content += '(\n'
            for ci in self.cell_info:
                self.formatted_content += "{} ".format(ci)
            self.formatted_content += '\n)'
        else:
            if self.cell_type == 0:
                self.formatted_content = "({} ({} {} {} {})".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], self.cell_type)  # declaration
            else:
                self.formatted_content = "({} ({} {} {} {} {})".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], self.cell_type, self.elem_type)

        self.formatted_content += ')'

    @classmethod
    def declaration(cls, num):
        return cls(0, 1, num, 0, 0)


class XF_MSH(object):
    def __init__(self):
        self.section_list = []

    def add_section(self, section: XF_Section):
        self.section_list.append(section)

    def save(self, fn):
        msh = open(fn, 'w')
        for i in range(len(self.section_list)):
            self.section_list[i].build_content()
            msh.write("{}\n".format(self.section_list[i].content))

        msh.close()

    @staticmethod
    def pnt_index(i: int, j: int, U: int):
        return i * U + j + 1

    @staticmethod
    def cell_index(i: int, j: int, d: int, U: int, V: int):
        if d == 1:
            return 0 if i == U - 1 or j == V - 1 else j * (U - 1) + i + 1
        elif d == 2:
            return 0 if i == 0 or j == V - 1 else j * (U - 1) + i
        elif d == 3:
            return 0 if i == 0 or j == 0 else (j - 1) * (U - 1) + i
        elif d == 4:
            return 0 if i == U - 1 or j == 0 else (j - 1) * (U - 1) + i + 1
        else:
            raise ValueError("Invalid direction!")

    @classmethod
    def build_from_2d(cls, grid):
        """
        从2维结构网格构建Fluent msh文件
        :param grid: 2D structural grid
        :return: XF_MSH object
        """

        U, V, Dim = grid.shape
        Dim = 2

        NodeCnt = U * V
        EdgeCnt = 2 * U * V - U - V
        CellCnt = (U - 1) * (V - 1)

        k = 0
        NodeList = np.empty((U * V, Dim), float)
        for i in range(U):
            for j in range(V):
                for d in range(Dim):
                    NodeList[k][d] = grid[i][j][d]
                k += 1

        msh = cls()
        msh.add_section(XF_Header())
        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(Dim))
        msh.add_section(XF_Comment("Declaration:"))
        msh.add_section(XF_Node.declaration(NodeCnt, Dim))
        msh.add_section(XF_Face.declaration(EdgeCnt))
        msh.add_section(XF_Cell.declaration(CellCnt))
        msh.add_section(XF_Comment("Grid:"))
        msh.add_section(XF_Node(1, 1, NodeCnt, 1, Dim, NodeList))
        msh.add_section(XF_Cell(2, 1, CellCnt, 1, 3))
        return msh
