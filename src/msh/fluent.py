import numpy as np
from abc import abstractmethod
from enum import Enum, unique


class XF_Section(object):
    def __init__(self, index: int):
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
    def __init__(self, info='Grid generated with Python 3.x'):
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


@unique
class NodeType(Enum):
    Virtual, Any, Boundary = range(3)


class XF_Node(XF_Section):
    def __init__(self, zone: int, first: int, last: int, tp: NodeType, dim: int, pts=None):
        """
        网格点
        :param zone: 区域号 
        :param first: 网格点起始序号(Starting from 1)
        :param last:  网格点终止序号
        :param tp: 网格点类型(TGrid usage)
        :param dim: Indicate the dimensionality of the node data.
        :param pts: 网格点数组，(last-first+1)个元素
        """

        super(XF_Node, self).__init__(10)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.node_type = tp.value
        self.ND = dim
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
            self.formatted_content += "("
            for i in range(n):
                for d in range(self.ND):
                    self.formatted_content += "{}{}".format('\n' if d == 0 else ' ', self.pts[i][d])
            self.formatted_content += ")"
        self.formatted_content += ")"

    @classmethod
    def declaration(cls, num, dim):
        return cls(0, 1, num, NodeType.Virtual, dim)


class BCType(Enum):
    Inactive = 0
    Interior = 2
    Wall = 3
    PressureInlet = 4
    InletVent = 4
    IntakeFan = 4
    PressureOutlet = 5
    OutletVent = 5
    ExhaustFan = 5
    Symmetry = 7
    PeriodicShadow = 8
    PressureFarField = 9
    VelocityInlet = 10
    Periodic = 12
    Fan = 14
    PorousJump = 14
    Radiator = 14
    MassFlowInlet = 20
    Interface = 24
    Parent = 31  # hanging node
    Outflow = 36
    Axis = 37


@unique
class FaceType(Enum):
    Mixed = 0
    Linear = 2
    Triangular = 3
    Quadrilateral = 4
    Polygonal = 5


class XF_Face(XF_Section):
    def __init__(self, zone: int, first: int, last: int, bct: BCType, ft: FaceType, face_info=None):
        """
        边界
        :param zone: 区域号
        :param first: 起始序号(Starting from 1)
        :param last: 终止序号
        :param bct: 边界属性
        :param ft: 形状类别
        :param face_info: 边界信息，一维数组(last-first+1)个元素，每个元素中包含了邻接关系
        """

        super(XF_Face, self).__init__(13)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.bc_type = bct.value
        self.face_type = ft.value
        self.face_info = None if face_info is None else np.copy(face_info)

    def build_content(self):
        if self.face_info is None:
            self.formatted_content = "({} ({} {} {} {})".format(self.index,
                                                                hex(self.zone_id)[2:],
                                                                hex(self.first_index)[2:],
                                                                hex(self.last_index)[2:],
                                                                hex(self.face_type)[2:])
        else:
            self.formatted_content = "({} ({} {} {} {} {})".format(self.index,
                                                                   hex(self.zone_id)[2:],
                                                                   hex(self.first_index)[2:],
                                                                   hex(self.last_index)[2:],
                                                                   hex(self.bc_type)[2:],
                                                                   hex(self.face_type)[2:])
            self.formatted_content += '('
            for fc in self.face_info:
                for cfi in range(len(fc)):  # current face info
                    self.formatted_content += "{}{}".format('\n' if cfi == 0 else ' ', fc[cfi])
            self.formatted_content += ')'

        self.formatted_content += ')'

    @classmethod
    def declaration(cls, num):
        return cls(0, 1, num, BCType.Inactive, FaceType.Mixed)


@unique
class CellType(Enum):
    Dead = 0
    Fluid = 1
    Solid = 17
    Inactive = 32


@unique
class CellElement(Enum):
    Mixed = 0
    Triangular = 1
    Tetrahedral = 2
    Quadrilateral = 3
    Hexahedra = 4
    Pyramid = 5
    Wedge = 6
    Polyhedral = 7


class XF_Cell(XF_Section):
    def __init__(self, zone: int, first: int, last: int, ct: CellType, ce: CellElement, cell_info=None):
        """
        单元
        :param zone: 区域号
        :param first: 起始序号(Starting from 1)
        :param last: 终止序号
        :param ct: Cell type.
        :param ce: Cell element.
        :param cell_info: 单元类型信息，一维数组，(last-first+1)个元素
        """

        super(XF_Cell, self).__init__(12)
        self.zone_id = zone
        self.first_index = first
        self.last_index = last
        self.cell_type = ct.value
        self.elem_type = ce.value
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
        return cls(0, 1, num, CellType.Dead, CellElement.Mixed)


class XF_MSH(object):
    def __init__(self):
        self.section_list = []
        self.content = ''

    def add_section(self, section: XF_Section):
        self.section_list.append(section)
        section.build_content()
        self.content += "{}\n".format(section.content)

    def add_blank(self):
        self.content += '\n'

    def save(self, fn):
        msh = open(fn, 'w')
        msh.write(self.content)
        msh.close()

    @staticmethod
    def pnt_index(i: int, j: int, U: int):
        return j * U + i + 1

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

    @staticmethod
    def intersect(p1, p2, U: int, V: int):
        i, j = p1
        ii, jj = p2
        di = ii - i
        dj = jj - j

        if di == 0:
            if dj == 1:
                return XF_MSH.cell_index(i, j, 2, U, V), XF_MSH.cell_index(i, j, 1, U, V)
            elif dj == -1:
                return XF_MSH.cell_index(i, j, 4, U, V), XF_MSH.cell_index(i, j, 3, U, V)
            else:
                raise ValueError("Invalid coordinates!")
        elif dj == 0:
            if di == 1:
                return XF_MSH.cell_index(i, j, 1, U, V), XF_MSH.cell_index(i, j, 4, U, V)
            elif di == -1:
                return XF_MSH.cell_index(i, j, 3, U, V), XF_MSH.cell_index(i, j, 2, U, V)
            else:
                raise ValueError("Invalid coordinates!")
        else:
            raise ValueError("Invalid coordinates!")

    @classmethod
    def from_str2d(cls, grid, bc=(BCType.Wall, BCType.VelocityInlet, BCType.Wall, BCType.Outflow)):
        """
        从2维结构网格构建Fluent msh文件
        :param grid: 2D structural grid
        :param bc: Boundary ID
        :return: XF_MSH object
        """

        msh = cls()
        U, V, Dim = grid.shape
        Dim = 2

        msh.add_section(XF_Header())
        msh.add_blank()
        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(Dim))
        msh.add_blank()

        NodeCnt = U * V
        EdgeCnt = 2 * U * V - U - V
        CellCnt = (U - 1) * (V - 1)

        msh.add_section(XF_Comment("Declaration:"))
        msh.add_section(XF_Cell.declaration(CellCnt))
        msh.add_section(XF_Face.declaration(EdgeCnt))
        msh.add_section(XF_Node.declaration(NodeCnt, Dim))
        msh.add_blank()

        '''Node'''
        k = 0
        NodeList = np.empty((NodeCnt, Dim), float)
        for i in range(U):
            for j in range(V):
                for d in range(Dim):
                    NodeList[k][d] = grid[i][j][d]
                k += 1

        msh.add_section(XF_Comment("Grid:"))
        msh.add_section(XF_Node(1, 1, NodeCnt, NodeType.Any, Dim, NodeList))
        msh.add_blank()

        '''Cell'''
        msh.add_section(XF_Cell(2, 1, CellCnt, CellType.Fluid, CellElement.Quadrilateral))
        msh.add_blank()

        '''Face-C1'''
        fcs = 1
        fc = U - 1
        Face1 = np.empty((fc, 4), int)
        for i in range(fc):
            Face1[i][0] = XF_MSH.pnt_index(i, 0, U)
            Face1[i][1] = XF_MSH.pnt_index(i + 1, 0, U)
            Face1[i][3], Face1[i][2] = XF_MSH.intersect((i, 0), (i + 1, 0), U, V)

        msh.add_section(XF_Face(3, fcs, fcs + fc - 1, bc[0], FaceType.Linear, Face1))
        msh.add_blank()

        '''Face-C2'''
        fcs = 1
        fc = V - 1
        Face2 = np.empty((fc, 4), int)
        for j in range(fc):
            Face2[j][0] = XF_MSH.pnt_index(0, j, U)
            Face2[j][1] = XF_MSH.pnt_index(0, j + 1, U)
            Face2[j][3], Face2[j][2] = XF_MSH.intersect((0, j), (0, j + 1), U, V)

        msh.add_section(XF_Face(4, fcs, fcs + fc - 1, bc[1], FaceType.Linear, Face2))
        msh.add_blank()

        '''Face-C3'''
        fcs = 1
        fc = U - 1
        Face3 = np.empty((fc, 4), int)
        for i in range(fc):
            Face3[i][0] = XF_MSH.pnt_index(i, V - 1, U)
            Face3[i][1] = XF_MSH.pnt_index(i + 1, V - 1, U)
            Face3[i][3], Face3[i][2] = XF_MSH.intersect((i, V - 1), (i + 1, V - 1), U, V)

        msh.add_section(XF_Face(5, fcs, fcs + fc - 1, bc[2], FaceType.Linear, Face3))
        msh.add_blank()

        '''Face-C4'''
        fcs = 1
        fc = V - 1
        Face4 = np.empty((fc, 4), int)
        for j in range(fc):
            Face4[j][0] = XF_MSH.pnt_index(U - 1, j, U)
            Face4[j][1] = XF_MSH.pnt_index(U - 1, j + 1, U)
            Face4[j][3], Face4[j][2] = XF_MSH.intersect((U - 1, j), (U - 1, j + 1), U, V)

        msh.add_section(XF_Face(6, fcs, fcs + fc - 1, bc[3], FaceType.Linear, Face4))
        msh.add_blank()

        '''Interior'''
        fcs = 1
        Face5 = []
        for i in range(1, U - 1):
            for j in range(V - 1):
                n1 = XF_MSH.pnt_index(i, j, U)
                n2 = XF_MSH.pnt_index(i, j + 1, U)
                cl, cr = XF_MSH.intersect((i, j), (i, j + 1), U, V)
                Face5.append(np.array([n1, n2, cr, cl], int))

        for j in range(1, V - 1):
            for i in range(U - 1):
                n1 = XF_MSH.pnt_index(i, j, U)
                n2 = XF_MSH.pnt_index(i + 1, j, U)
                cl, cr = XF_MSH.intersect((i, j), (i + 1, j), U, V)
                Face5.append(np.array([n1, n2, cr, cl], int))

        fc = len(Face5)
        Face5 = np.copy(Face5)

        msh.add_section(XF_Face(7, fcs, fc, BCType.Interior, FaceType.Linear, Face5))

        return msh
