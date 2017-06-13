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
        self.formatted_content = "({} ({} {} {} {} {})".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.node_type)[2:], self.ND)
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
        if self.face_info is not None:
            self.formatted_content = "({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.bc_type)[2:], hex(self.face_type)[2:])
            for fc in self.face_info:
                for cfi in range(len(fc)):  # current face info
                    self.formatted_content += "{}{}".format('\n' if cfi == 0 else ' ', hex(fc[cfi])[2:])
            self.formatted_content += '))'
        else:
            self.formatted_content = "({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.face_type)[2:])

    @classmethod
    def declaration(cls, num):
        return cls(0, 1, num, BCType.Inactive, FaceType.Mixed)


@unique
class CellType(Enum):
    Dead = 0
    Fluid = 1
    Solid = 17


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
            self.formatted_content = "({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:], hex(self.elem_type)[2:])
            self.formatted_content += "\n{}".format(hex(self.cell_info[0])[2:])
            for ci in range(1, len(self.cell_info)):
                self.formatted_content += " {}".format(hex(self.cell_info[ci])[2:])
            self.formatted_content += "))"
        else:
            if self.cell_type == 0:
                self.formatted_content = "({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:])  # declaration
            else:
                self.formatted_content = "({} ({} {} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:], hex(self.elem_type)[2:])

    @classmethod
    def declaration(cls, num):
        return cls(0, 1, num, CellType.Dead, CellElement.Mixed)


class XF_MSH(object):
    def __init__(self):
        """
        ANSYS Fluent MSH文件
        """

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
        """
        坐标对应的序号，序号从1开始
        :param i: 列号(Starting from 0)，与x方向对应
        :param j: 行号(Starting from 0)，与y方向对应
        :param U: 每行节点数
        :return: (i,j)节点所对应的序号 
        """

        return j * U + i + 1

    @staticmethod
    def cell_index(i: int, j: int, d: int, U: int, V: int):
        """
        由节点坐标及象限编号确定单元序号
        :param i: 列号(Starting from 0)，与x方向对应
        :param j: 行号(Starting from 0)，与y方向对应
        :param d: 象限号(Range from 1 to 4)
        :param U: 每行节点数
        :param V: 每列节点数
        :return: 以(i,j)为原点，位于第d象限的Cell的序号(Starting from 1)，若不存在则为0
        """

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
        """
        求两个相邻的节点连结成的线段左右两个Cell
        :param p1: 起始点坐标
        :param p2: 终止点坐标
        :param U: 每行节点数
        :param V: 每列节点数
        :return: 从p1到p2的向量左右两个Cell的序号: cl, cr
        """

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

        U, V, Dim = grid.shape
        Dim = 2
        NodeCnt = U * V
        EdgeCnt = 2 * U * V - U - V
        CellCnt = (U - 1) * (V - 1)

        k = 0
        NodeList = np.empty((NodeCnt, Dim), float)
        for j in range(V):
            for i in range(U):
                for d in range(Dim):
                    NodeList[k][d] = grid[i][j][d]  # 节点序号规定为沿x方向优先，依次递增
                k += 1

        '''Face-C1'''
        face1_num = U - 1
        face1 = np.empty((face1_num, 4), int)
        face1_first = 1
        face1_last = face1_first + face1_num - 1
        for i in range(face1_num):
            face1[i][0] = XF_MSH.pnt_index(i, 0, U)
            face1[i][1] = XF_MSH.pnt_index(i + 1, 0, U)
            cl, cr = XF_MSH.intersect((i, 0), (i + 1, 0), U, V)
            face1[i][2] = cl
            face1[i][3] = cr

        '''Face-C2'''
        face2_num = V - 1
        face2 = np.empty((face2_num, 4), int)
        face2_first = face1_last + 1
        face2_last = face2_first + face2_num - 1
        for j in range(face2_num):
            face2[j][0] = XF_MSH.pnt_index(0, j, U)
            face2[j][1] = XF_MSH.pnt_index(0, j + 1, U)
            cl, cr = XF_MSH.intersect((0, j), (0, j + 1), U, V)
            face2[j][2] = cl
            face2[j][3] = cr

        '''Face-C3'''
        face3_num = U - 1
        face3 = np.empty((face3_num, 4), int)
        face3_first = face2_last + 1
        face3_last = face3_first + face3_num - 1
        for i in range(face3_num):
            face3[i][0] = XF_MSH.pnt_index(i, V - 1, U)
            face3[i][1] = XF_MSH.pnt_index(i + 1, V - 1, U)
            cl, cr = XF_MSH.intersect((i, V - 1), (i + 1, V - 1), U, V)
            face3[i][2] = cl
            face3[i][3] = cr

        '''Face-C4'''
        face4_num = V - 1
        face4 = np.empty((face4_num, 4), int)
        face4_first = face3_last + 1
        face4_last = face4_first + face4_num - 1
        for j in range(face4_num):
            face4[j][0] = XF_MSH.pnt_index(U - 1, j, U)
            face4[j][1] = XF_MSH.pnt_index(U - 1, j + 1, U)
            cl, cr = XF_MSH.intersect((U - 1, j), (U - 1, j + 1), U, V)
            face4[j][2] = cl
            face4[j][3] = cr

        '''Interior'''
        face5 = []
        for i in range(1, U - 1):
            for j in range(V - 1):
                n1 = XF_MSH.pnt_index(i, j, U)
                n2 = XF_MSH.pnt_index(i, j + 1, U)
                cl, cr = XF_MSH.intersect((i, j), (i, j + 1), U, V)
                face5.append(np.array([n1, n2, cl, cr], int))

        for j in range(1, V - 1):
            for i in range(U - 1):
                n1 = XF_MSH.pnt_index(i, j, U)
                n2 = XF_MSH.pnt_index(i + 1, j, U)
                cl, cr = XF_MSH.intersect((i, j), (i + 1, j), U, V)
                face5.append(np.array([n1, n2, cl, cr], int))

        face5_num = len(face5)
        face5 = np.copy(face5)
        face5_first = face4_last + 1
        face5_last = face5_first + face5_num - 1

        '''MSH File'''
        msh = cls()
        msh.add_section(XF_Header())
        msh.add_blank()
        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(Dim))
        msh.add_blank()
        msh.add_section(XF_Comment("Declaration:"))
        msh.add_section(XF_Cell.declaration(CellCnt))
        msh.add_section(XF_Face.declaration(EdgeCnt))
        msh.add_section(XF_Node.declaration(NodeCnt, Dim))
        msh.add_blank()
        msh.add_section(XF_Comment("Grid:"))
        msh.add_section(XF_Cell(7, 1, CellCnt, CellType.Fluid, CellElement.Quadrilateral))
        msh.add_blank()
        msh.add_section(XF_Face(2, face1_first, face1_last, bc[0], FaceType.Linear, face1))
        msh.add_blank()
        msh.add_section(XF_Face(3, face2_first, face2_last, bc[1], FaceType.Linear, face2))
        msh.add_blank()
        msh.add_section(XF_Face(4, face3_first, face3_last, bc[2], FaceType.Linear, face3))
        msh.add_blank()
        msh.add_section(XF_Face(5, face4_first, face4_last, bc[3], FaceType.Linear, face4))
        msh.add_blank()
        msh.add_section(XF_Face(6, face5_first, face5_last, BCType.Interior, FaceType.Linear, face5))
        msh.add_blank()
        msh.add_section(XF_Node(1, 1, NodeCnt, NodeType.Any, Dim, NodeList))

        return msh
