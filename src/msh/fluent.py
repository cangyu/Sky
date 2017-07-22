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
        if self.pts is None:
            self.formatted_content = "({} ({} {} {} {})".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.node_type)[2:])
        else:
            self.formatted_content = "({} ({} {} {} {} {})".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.node_type)[2:], self.ND)
            n, dim = self.pts.shape
            self.formatted_content += "("
            for i in range(n):
                for d in range(self.ND):
                    self.formatted_content += "{}{}".format('\n' if d == 0 else ' ', self.pts[i][d])
            self.formatted_content += ")"
        self.formatted_content += ")"

    @classmethod
    def declaration(cls, num):
        return cls(0, 1, num, NodeType.Virtual, 0)


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
    Hexahedral = 4
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
    def dimensional_copy(dst, src, dim):
        for i in range(dim):
            dst[i] = src[i]

    @staticmethod
    def pnt_idx_2d(i, j, u):
        """
        坐标对应的序号，序号从1开始
        :param i: 列号(Starting from 0)，与x方向对应
        :type i: int
        :param j: 行号(Starting from 0)，与y方向对应
        :type j: int
        :param u: 每行节点数
        :type u: int
        :return: (i,j)节点所对应的序号 
        :rtype: int
        """

        return 1 + j * u + i

    @staticmethod
    def cell_idx_2d(i, j, d, u, v):
        """
        由节点坐标及象限编号确定单元序号
        :param i: 列号(Starting from 0)，与x方向对应
        :type i: int
        :param j: 行号(Starting from 0)，与y方向对应
        :type j: int
        :param d: 象限号(Range from 1 to 4)
        :type d: int
        :param u: 每行节点数
        :type u: int
        :param v: 每列节点数
        :type v: int
        :return: 以(i,j)为原点，位于第d象限的Cell的序号(Starting from 1)，若不存在则为0
        :rtype: int
        """

        if d == 1:
            return 0 if i == u - 1 or j == v - 1 else j * (u - 1) + i + 1
        elif d == 2:
            return 0 if i == 0 or j == v - 1 else j * (u - 1) + i
        elif d == 3:
            return 0 if i == 0 or j == 0 else (j - 1) * (u - 1) + i
        elif d == 4:
            return 0 if i == u - 1 or j == 0 else (j - 1) * (u - 1) + i + 1
        else:
            raise ValueError("Invalid direction!")

    @staticmethod
    def intersect(p1, p2, u, v):
        """
        求两个相邻的节点连结成的线段左右两个Cell
        :param p1: 起始点坐标
        :param p2: 终止点坐标
        :param u: 每行节点数
        :type u: int
        :param v: 每列节点数
        :type v: int
        :return: 从p1到p2的向量左右两个Cell的序号: cl, cr
        :rtype: (int, int)
        """

        i, j = p1
        ii, jj = p2
        di = ii - i
        dj = jj - j

        if di != 0 and dj != 0:
            raise AssertionError("Invalid coordinates!")

        if di == 0:
            if dj == 1:
                return XF_MSH.cell_idx_2d(i, j, 2, u, v), XF_MSH.cell_idx_2d(i, j, 1, u, v)
            elif dj == -1:
                return XF_MSH.cell_idx_2d(i, j, 4, u, v), XF_MSH.cell_idx_2d(i, j, 3, u, v)
            else:
                raise ValueError("Invalid coordinates!")
        else:
            if di == 1:
                return XF_MSH.cell_idx_2d(i, j, 1, u, v), XF_MSH.cell_idx_2d(i, j, 4, u, v)
            elif di == -1:
                return XF_MSH.cell_idx_2d(i, j, 3, u, v), XF_MSH.cell_idx_2d(i, j, 2, u, v)
            else:
                raise ValueError("Invalid coordinates!")

    @classmethod
    def from_str2d(cls, grid, bc=(BCType.Wall, BCType.VelocityInlet, BCType.Wall, BCType.Outflow)):
        """
        从2维结构网格构建Fluent msh文件
        :param grid: 2D structural grid
        :param bc: Boundary ID
        :return: 可用于后续生成msh文件的XF_MSH对象
        :rtype: XF_MSH
        """

        U, V, Dim = grid.shape
        Dim = 2
        NodeCnt = U * V
        EdgeCnt = 2 * U * V - U - V
        CellCnt = (U - 1) * (V - 1)
        zone_idx = 0

        k = 0
        NodeList = np.empty((NodeCnt, Dim), float)
        for j in range(V):
            for i in range(U):
                cls.dimensional_copy(NodeList[k], grid[i][j], Dim)  # 节点序号规定为沿x方向优先，依次递增
                k += 1

        '''Face-C1'''
        face1_num = U - 1
        face1 = np.empty((face1_num, 4), int)
        face1_first = 1
        face1_last = face1_first + face1_num - 1
        for i in range(face1_num):
            face1[i][0] = XF_MSH.pnt_idx_2d(i, 0, U)
            face1[i][1] = XF_MSH.pnt_idx_2d(i + 1, 0, U)
            cl, cr = XF_MSH.intersect((i, 0), (i + 1, 0), U, V)
            face1[i][2] = cl
            face1[i][3] = cr

        '''Face-C2'''
        face2_num = V - 1
        face2 = np.empty((face2_num, 4), int)
        face2_first = face1_last + 1
        face2_last = face2_first + face2_num - 1
        for j in range(face2_num):
            face2[j][0] = XF_MSH.pnt_idx_2d(0, j, U)
            face2[j][1] = XF_MSH.pnt_idx_2d(0, j + 1, U)
            cl, cr = XF_MSH.intersect((0, j), (0, j + 1), U, V)
            face2[j][2] = cl
            face2[j][3] = cr

        '''Face-C3'''
        face3_num = U - 1
        face3 = np.empty((face3_num, 4), int)
        face3_first = face2_last + 1
        face3_last = face3_first + face3_num - 1
        for i in range(face3_num):
            face3[i][0] = XF_MSH.pnt_idx_2d(i, V - 1, U)
            face3[i][1] = XF_MSH.pnt_idx_2d(i + 1, V - 1, U)
            cl, cr = XF_MSH.intersect((i, V - 1), (i + 1, V - 1), U, V)
            face3[i][2] = cl
            face3[i][3] = cr

        '''Face-C4'''
        face4_num = V - 1
        face4 = np.empty((face4_num, 4), int)
        face4_first = face3_last + 1
        face4_last = face4_first + face4_num - 1
        for j in range(face4_num):
            face4[j][0] = XF_MSH.pnt_idx_2d(U - 1, j, U)
            face4[j][1] = XF_MSH.pnt_idx_2d(U - 1, j + 1, U)
            cl, cr = XF_MSH.intersect((U - 1, j), (U - 1, j + 1), U, V)
            face4[j][2] = cl
            face4[j][3] = cr

        '''Interior'''
        face5 = []
        for i in range(1, U - 1):
            for j in range(V - 1):
                n1 = XF_MSH.pnt_idx_2d(i, j, U)
                n2 = XF_MSH.pnt_idx_2d(i, j + 1, U)
                cl, cr = XF_MSH.intersect((i, j), (i, j + 1), U, V)
                face5.append(np.array([n1, n2, cl, cr], int))

        for j in range(1, V - 1):
            for i in range(U - 1):
                n1 = XF_MSH.pnt_idx_2d(i, j, U)
                n2 = XF_MSH.pnt_idx_2d(i + 1, j, U)
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
        msh.add_section(XF_Node.declaration(NodeCnt))
        msh.add_blank()

        msh.add_section(XF_Comment("Grid:"))

        zone_idx += 1
        msh.add_section(XF_Node(zone_idx, 1, NodeCnt, NodeType.Any, Dim, NodeList))
        msh.add_blank()

        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, face1_first, face1_last, bc[0], FaceType.Linear, face1))
        msh.add_blank()

        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, face2_first, face2_last, bc[1], FaceType.Linear, face2))
        msh.add_blank()

        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, face3_first, face3_last, bc[2], FaceType.Linear, face3))
        msh.add_blank()

        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, face4_first, face4_last, bc[3], FaceType.Linear, face4))
        msh.add_blank()

        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, face5_first, face5_last, BCType.Interior, FaceType.Linear, face5))
        msh.add_blank()

        zone_idx += 1
        msh.add_section(XF_Cell(zone_idx, 1, CellCnt, CellType.Fluid, CellElement.Quadrilateral))

        return msh

    @classmethod
    def from_str2d_multi(cls, blk_list, bc_list, adj_info):
        """
        根据给定的邻接关系，将多块结构网格转换成非结构形式，
        并赋予给定的边界条件, 类似于ICEM中的'Convert to unstructured mesh'功能
        :param blk_list: 单块结构网格序列
        :param bc_list: 依次包括每个单块结构网格的BC
        :param adj_info: 依次包括每个单块结构网格的邻接关系序列
        :return: 可用于后续生成msh文件的XF_MSH对象
        :rtype: XF_MSH
        """

        def cell_num(u, v):
            """
            计算单个Block中的Cell数量
            :param u: U方向节点数
            :type u: int
            :param v: V方向节点数
            :type v: int
            :return: (U,V)方向上长度分别为(u,v)的Block的Cell数量
            :rtype: int
            """

            return (u - 1) * (v - 1)

        def edge_num(u, v):
            """
            计算单个Block中的Edge数量
            :param u: U方向节点数
            :type u: int
            :param v: V方向节点数
            :type v: int
            :return: (U,V)方向上长度分别为(u,v)的Block的Edge数量
            :rtype: int
            """

            return (u - 1) * v + (v - 1) * u

        def internal_edge_num(u, v):
            """
            Block内部Edge数量
            :param u: U方向节点数
            :type u: int
            :param v: V方向节点数
            :type v: int
            :return: (U,V)方向上长度分别为(u,v)的Block的内部所有Edge数量
            :rtype: int
            """

            return (u - 2) * (v - 1) + (v - 2) * (u - 1)

        def boundary_edge_num(u, v, e=None):
            """
            计算Block指定边上的Edge数量
            :param u: U方向节点数
            :type u: int
            :param v: V方向节点数
            :type v: int
            :param e: 边序号，从1开始，若为None则表示计算4条边界上的Edge数量
            :return: (U,V)方向上长度分别为(u,v)的Block的第e条边上的Edge数量
            """

            if e is None:
                return 2 * (u + v - 2)  # all boundary edge num
            elif e in (1, 3):
                return u - 1
            elif e in (2, 4):
                return v - 1
            else:
                raise ValueError("Invalid edge index {}".format(e))

        def blk_edge_pnt(blk_idx, edge):
            """
            给定Block序号和边的序号，计算其上节点数量
            :param blk_idx: Block序号，从0开始
            :type blk_idx: int
            :param edge: 边的序号，从1开始
            :type edge: int
            :return: 第blk_idx个Block的第edge条边上节点的数量
            :rtype: int
            """

            return blk_shape[blk_idx][0] if edge in (1, 3) else blk_shape[blk_idx][1]

        def is_boundary_pnt(i, j, u, v):
            """
            判断点是否是边界点
            :param i: 目标点U方向下标
            :type i: int
            :param j: 目标点V方向下标
            :type j: int
            :param u: U方向节点数
            :type u: int
            :param v: V方向节点数
            :type v: int
            :return: 在(U,V)方向长度分别为(u,v)的Block的边界上，点(i,j)是否在其边界上
            :rtype: bool
            """

            return i in (0, u - 1) or j in (0, v - 1)

        def boundary_pnt_num(u, v):
            """
            计算Block最外层上的节点数量
            :param u: U方向节点数
            :type u: int
            :param v: V方向节点数
            :type v: int
            :return: Block最外层上的节点数量
            :rtype: int
            """

            return 2 * (u + v - 2)

        def boundary_pnt_idx(i, j, u, v):
            """
            根据坐标计算边界点的序号
            :param i: 目标点U方向下标
            :type i: int
            :param j: 目标点V方向下标
            :type j: int
            :param u: U方向节点数
            :type u: int
            :param v: V方向节点数
            :type v: int
            :return: 在(U,V)方向长度分别为(u,v)的Block的边界上，点(i,j)的序号，从0开始
            :rtype: int
            """

            if i == 0 and j == 0:
                return 0
            elif i == u - 1 and j == 0:
                return u - 1
            elif i == u - 1 and j == v - 1:
                return u + v - 2
            elif i == 0 and j == v - 1:
                return 2 * u + v - 3
            else:
                idx = get_edge_idx(i, j, u, v)[0]
                if idx == 1:
                    return i
                elif idx == 4:
                    return u - 1 + j
                elif idx == 3:
                    return 2 * u + v - 3 - i
                elif idx == 2:
                    return 2 * (u + v) - 4 - j
                else:
                    raise ValueError('Invalid edge index: {}'.format(idx))

        def get_pnt_idx(i, j, k):
            """
            求给定点在本程序的编号规则下的序号
            :param i: 目标点U方向下标
            :type i: int
            :param j: 目标点V方向下标
            :type j: int
            :param k: Block序号
            :type k: int
            :return: 第(k+1)个Block中点(i,j)在本程序的编号规则下的序号
            :rtype: int
            """

            u, v = blk_shape[k]
            if is_boundary_pnt(i, j, u, v):
                idx = boundary_pnt_idx(i, j, u, v)
                return bc_flag[k][idx]
            else:
                rel = (j - 1) * (u - 2) + (i - 1)
                return internal_pnt_start[k] + rel

        def boundary_coordinate(idx, u, v):
            """
            给定边界点序号反算其坐标
            边界点序号定义为从原点开始，按逆时针方向依次遍历4条边界上的所有点
            :param idx: 边界点序号，从0开始
            :type idx: int
            :param u: U方向节点数
            :type u: int
            :param v: V方向节点数
            :type v: int
            :return: 在(U,V)方向上长度分别为(u,v)的Block上，其边界上第(idx+1)个点的坐标
            """

            s1, s2, s3, s4 = u - 1, u + v - 2, 2 * u + v - 3, boundary_pnt_num(u, v)

            if idx < 0 or idx >= s4:
                raise ValueError("Invalid boundary point index: {}".format(idx))

            if idx < s1:
                return idx, 0
            elif idx < s2:
                return s1, idx - s1
            elif idx < s3:
                return s3 - idx, v - 1
            else:
                return 0, s4 - idx

        def get_edge_idx(i, j, u, v):
            """
            根据坐标计算目标点在哪条边上
            若是角点则有2条边，若是其它边界点则有1条边，否则为空
            :param i: 目标点U方向下标
            :type i: int
            :param j: 目标点V方向下标
            :type j: int
            :param u: U方向节点数
            :type u: int
            :param v: V方向节点数
            :type v: int
            :return: 在(U,V)方向上长度分别为(u,v)的Block上，点(i,j)所处边的序号列表，从1开始
            """

            ans = []

            if i == 0:
                ans.append(2)
            if i == u - 1:
                ans.append(4)
            if j == 0:
                ans.append(1)
            if j == v - 1:
                ans.append(3)

            return ans

        def get_counterpart_idx(k1, e1, i, j, k2, e2):
            """
            给定Block1上的点和该点对应的边序号，根据邻接关系求出与该边对应的Block2上对应边上的点在Block2边界上的序号
            先计算对应点的坐标，再将坐标转换成序号
            :param k1: Block1序号，从0开始
            :type k1: int
            :param e1: Block1上边界点所在的边序号，从1开始
            :type k1: int
            :param i: Block1上边界点U方向下标
            :type k1: int
            :param j: Block1上边界点V方向下标
            :type k1: int
            :param k2: Block1序号，从0开始
            :type k1: int
            :param e2: Block2上对应边的序号，从1开始
            :type k1: int
            :return: 给定第(k1+1)个Block上第e1条边上的点(i,j), 计算在第(k2+1)个Block的第e2条边上所对应的点在该Block上的序号
            :rtype: int
            """

            u2, v2 = blk_shape[k2]

            if e1 in (1, 3):
                t = i
            elif e1 in (2, 4):
                t = j
            else:
                raise ValueError('Invalid edge index: {}'.format(e1))

            if e2 == 1:
                ii, jj = t, 0
            elif e2 == 2:
                ii, jj = 0, t
            elif e2 == 3:
                ii, jj = t, v2 - 1
            elif e2 == 4:
                ii, jj = u2 - 1, t
            else:
                raise ValueError("Invalid counterpart edge index: {}".format(e2))

            return boundary_pnt_idx(ii, jj, u2, v2)

        def boundary_cell_list(k, e):
            """
            某条边界上的Cell列表
            :param k: Block序号，从0开始
            :type k: int
            :param e: 边序号，从1开始
            :type e: int
            :return: 第(k+1)个Block的第e条边上的Cell序号
            """

            u, v = blk_shape[k]
            if e == 1:
                start = cell_start[k]
                end = start + (u - 1)
                gap = 1
            elif e == 2:
                start = cell_start[k]
                end = start + (u - 1) * (v - 1)
                gap = u - 1
            elif e == 3:
                start = cell_start[k] + (u - 1) * (v - 2)
                end = start + (u - 1)
                gap = 1
            elif e == 4:
                start = cell_start[k] + (u - 2)
                end = cell_start[k] + (u - 1) * (v - 1)
                gap = u - 1
            else:
                raise ValueError("Invalid edge index!")

            return np.arange(start, end, gap, dtype=int)

        def boundary_pnt_idx_list(k, e):
            """
            某条边界上节点的序号列表
            :param k: Block序号，从0开始
            :type k: int
            :param e: 边序号，从1开始
            :type e: int
            :return: 第(k+1)个Block的第e条边上所有节点的序号
            """

            u, v = blk_shape[k]
            if e == 1:
                return bc_flag[k][:u]
            elif e == 2:
                ans = np.array([bc_flag[k][0]])
                tmp = bc_flag[k][2 * u + v - 3:]
                tmp = tmp[::-1]
                ans = np.append(ans, tmp)
                return ans
            elif e == 3:
                ans = bc_flag[k][u + v - 2:2 * u + v - 2]
                ans = ans[::-1]
                return ans
            elif e == 4:
                return bc_flag[k][u - 1:u + v - 1]
            else:
                raise ValueError("Invalid edge index!")

        def handle_interior_edge(k1, e1, k2, e2):
            """
            构建重合的Interior边的描述数组
            :param k1: Block1的序号，从0开始
            :type k1: int
            :param e1: 该Interior边在Block1中的边序号，从1开始
            :type e1: int
            :param k2: Block2的序号，从0开始
            :type k2: int
            :param e2: 该Interior边在Block2中的边序号，从1开始
            :type e2: int
            :return: 表示该Interior边上所有Edge的邻接关系的数组
            """

            bc1 = boundary_cell_list(k1, e1)
            bc2 = boundary_cell_list(k2, e2)

            cn = blk_edge_pnt(k1, e1) - 1
            cur_boundary_edge = np.empty((cn, 4), int)
            pnt_idx_list = boundary_pnt_idx_list(k1, e1)

            for i in range(cn):
                n1 = pnt_idx_list[i]
                n2 = pnt_idx_list[i + 1]
                cl = bc1[i]
                cr = bc2[i]
                cur_boundary_edge[i] = np.array([n1, n2, cl, cr], int)

            return cur_boundary_edge

        @unique
        class DirRevIdFactor(Enum):
            Edge1Reverted = 2
            Edge2Reverted = 3
            Edge3Reverted = 5
            Edge4Reverted = 7
            EdgeRevTester = 210

        '''Basic variables'''
        dimension = 2
        total_blk = len(blk_list)
        blk_shape = np.empty((total_blk, 2), int)
        adj_desc = np.zeros((total_blk, 4, 2), int)
        blk_rev = np.ones(total_blk, int)
        for k, blk in enumerate(blk_list):
            blk_shape[k] = blk.shape[:2]

        '''Adjacent info and orientation'''
        for entry in adj_info:
            b1, e1 = entry[0]
            b2, e2 = entry[1]
            if b1 >= total_blk or b2 >= total_blk or e1 > 4 or e2 > 4:
                raise ValueError("Invalid input.")

            if e1 != 0 and e2 != 0:
                adj_desc[b1][e1 - 1] = np.array([b2, e2], int)
                adj_desc[b2][e2 - 1] = np.array([b1, e1], int)

            if e1 == 2:
                blk_rev[b1] *= DirRevIdFactor.Edge2Reverted.value
            elif e1 == 3:
                blk_rev[b1] *= DirRevIdFactor.Edge3Reverted.value
            else:
                pass

            if e2 == 1:
                blk_rev[b2] *= DirRevIdFactor.Edge1Reverted.value
            elif e2 == 4:
                blk_rev[b2] *= DirRevIdFactor.Edge4Reverted.value
            else:
                pass

        '''Define orientation'''
        for i in range(total_blk):
            if blk_rev[i] == 1:
                blk_rev[i] = 0
            else:
                blk_rev[i] = 1

        '''Initialize MSH file'''
        msh = cls()
        zone_idx = 0
        msh.add_section(XF_Header())
        msh.add_blank()

        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(dimension))
        msh.add_blank()

        '''Counting cells'''
        cell_start = np.empty(total_blk + 1, int)
        cell_start[0] = 1
        for k in range(total_blk):
            cu, cv = blk_shape[k]
            cur_cell_cnt = cell_num(cu, cv)
            cell_start[k + 1] = cell_start[k] + cur_cell_cnt
        total_cell = cell_start[-1] - 1

        '''Flush cell info to MSH file'''
        msh.add_section(XF_Comment("Cell Declaration:"))
        msh.add_section(XF_Cell.declaration(total_cell))
        msh.add_blank()

        zone_idx += 1
        msh.add_section(XF_Comment("Cell Info:"))
        msh.add_section(XF_Cell(zone_idx, 1, total_cell, CellType.Fluid, CellElement.Quadrilateral))
        msh.add_blank()

        '''Counting points'''
        total_pnt = 0
        for k in range(total_blk):
            cu, cv = blk_shape[k]
            total_pnt += cu * cv

        for entry in adj_info:
            b1, e1 = entry[0]
            b2, e2 = entry[1]
            if e1 != 0 and e2 != 0:
                cur_pnt_num = blk_edge_pnt(b1, e1)
                total_pnt -= cur_pnt_num

        '''Flush point declaration msg to MSH file'''
        msh.add_section(XF_Comment("Point Declaration:"))
        msh.add_section(XF_Node.declaration(total_pnt))
        msh.add_blank()

        '''Local storage for recording boundary point index'''
        bc_flag = []
        for k, blk in enumerate(blk_list):
            cu, cv = blk_shape[k]
            cur_pnt_num = boundary_pnt_num(cu, cv)
            flag_ary = np.full(cur_pnt_num, -1, int)  # -1表示该点未分配序号
            bc_flag.append(flag_ary)

        pnt_list = np.empty((total_pnt, dimension), float)
        pnt_idx = 0

        '''Handling internal points first'''
        internal_pnt_start = np.empty(total_blk + 1, int)
        internal_pnt_start[0] = 1
        for k, blk in enumerate(blk_list):
            cu, cv = blk_shape[k]
            for j in range(1, cv - 1):
                for i in range(1, cu - 1):
                    cls.dimensional_copy(pnt_list[pnt_idx], blk[i][j], dimension)
                    pnt_idx += 1
            internal_pnt_start[k + 1] = pnt_idx + 1

        '''Handling boundary points afterwards'''
        for k, blk in enumerate(blk_list):
            cu, cv = blk_shape[k]
            cn = boundary_pnt_num(cu, cv)
            for w in range(cn):
                if bc_flag[k][w] != -1:  # Has been marked
                    continue

                i, j = boundary_coordinate(w, cu, cv)
                edge_idx = get_edge_idx(i, j, cu, cv)  # Usually 1 edge, but there will be 2 edges when it's a corner point

                '''Check if it's an interior pnt'''
                is_interior = False
                interior_adj_edge = []
                for e in edge_idx:
                    if adj_desc[k][e - 1][1] != 0:  # 非0则表示有相邻的边
                        interior_adj_edge.append(e)
                        if not is_interior:
                            is_interior = True

                if not is_interior:
                    cls.dimensional_copy(pnt_list[pnt_idx], blk[i][j], dimension)
                    pnt_idx += 1
                    bc_flag[k][w] = pnt_idx
                else:
                    has_assigned = False
                    cur_adj = []

                    '''Inspect neighbours to see if it has been marked or not'''
                    for e in interior_adj_edge:
                        k2, e2 = adj_desc[k][e - 1]
                        cp_idx = get_counterpart_idx(k, e, i, j, k2, e2)
                        cur_adj.append((k2, cp_idx))
                        if bc_flag[k2][cp_idx] != -1:  # Neighbours has been marked
                            bc_flag[k][w] = bc_flag[k2][cp_idx]
                            has_assigned = True
                            break

                    if not has_assigned:
                        cls.dimensional_copy(pnt_list[pnt_idx], blk[i][j], dimension)
                        pnt_idx += 1
                        bc_flag[k][w] = pnt_idx
                        for e in cur_adj:
                            bc_flag[e[0]][e[1]] = pnt_idx

        '''Flush point coordinates to MSH file'''
        zone_idx += 1
        msh.add_section(XF_Comment("Point Coordinates:"))
        msh.add_section(XF_Node(zone_idx, 1, total_pnt, NodeType.Any, dimension, pnt_list))
        msh.add_blank()

        '''Counting edges'''
        total_edge = 0
        for k in range(total_blk):
            cu, cv = blk_shape[k]
            total_edge += edge_num(cu, cv)

        for entry in adj_info:
            k1, e1 = entry[0]
            k2, e2 = entry[1]
            if e1 != 0 and e2 != 0:
                cu, cv = blk_shape[k1]
                total_edge -= boundary_edge_num(cu, cv, e1)

        '''Flush edge declaration to MSH file'''
        msh.add_section(XF_Comment("Edge Declaration:"))
        msh.add_section(XF_Face.declaration(total_edge))
        msh.add_blank()

        '''Handle internal edges in each blk'''
        cur_edge_first = 1
        for k, blk in enumerate(blk_list):
            cu, cv = blk_shape[k]
            cn = internal_edge_num(cu, cv)
            cur_inter_edge = np.empty((cn, 4), int)
            ce = 0
            ts = cell_start[k] - 1

            '''Horizontal'''
            for j in range(1, cv - 1):
                for i in range(cu - 1):
                    n1 = get_pnt_idx(i, j, k)
                    n2 = get_pnt_idx(i + 1, j, k)
                    cl, cr = XF_MSH.intersect((i, j), (i + 1, j), cu, cv)
                    cl += ts
                    cr += ts
                    if blk_rev[k]:
                        cl, cr = cr, cl
                    cur_inter_edge[ce] = np.array([n1, n2, cl, cr], int)
                    ce += 1

            '''Vertical'''
            for i in range(1, cu - 1):
                for j in range(cv - 1):
                    n1 = get_pnt_idx(i, j, k)
                    n2 = get_pnt_idx(i, j + 1, k)
                    cl, cr = XF_MSH.intersect((i, j), (i, j + 1), cu, cv)
                    cl += ts
                    cr += ts
                    if blk_rev[k]:
                        cl, cr = cr, cl
                    cur_inter_edge[ce] = np.array([n1, n2, cl, cr], int)
                    ce += 1

            '''Flush internal edges into MSH file'''
            next_edge_first = cur_edge_first + ce
            msh.add_section(XF_Comment("Block {} internal edges:".format(k)))
            zone_idx += 1
            msh.add_section(XF_Face(zone_idx, cur_edge_first, next_edge_first - 1, BCType.Interior, FaceType.Linear, cur_inter_edge))
            msh.add_blank()
            cur_edge_first = next_edge_first

        '''Handle interior boundary'''
        for k, entry in enumerate(adj_info):
            k1, e1 = entry[0]
            k2, e2 = entry[1]
            if e1 == 0 or e2 == 0:
                continue

            cur_edge_list = handle_interior_edge(k1, e1, k2, e2)

            '''Flush interior edges into MSH file'''
            next_edge_first = cur_edge_first + len(cur_edge_list)
            msh.add_section(XF_Comment("Interior edge {}:".format(k)))
            zone_idx += 1
            msh.add_section(XF_Face(zone_idx, cur_edge_first, next_edge_first - 1, BCType.Interior, FaceType.Linear, cur_edge_list))
            msh.add_blank()
            cur_edge_first = next_edge_first

        '''Handle non-adjacent boundary'''
        be = 0
        for k in range(total_blk):
            for e in range(1, 5):
                if adj_desc[k][e - 1][1] == 0:
                    be += 1
                    cur_edge_cell = boundary_cell_list(k, e)
                    cur_edge_pnt = boundary_pnt_idx_list(k, e)
                    cn = blk_edge_pnt(k, e) - 1
                    cur_edge_list = np.empty((cn, 4), int)

                    for i in range(cn):
                        n1 = cur_edge_pnt[i]
                        n2 = cur_edge_pnt[i + 1]
                        cl = cur_edge_cell[i]
                        cr = 0
                        if (blk_rev[k] == 0 and e in (2, 3)) or (blk_rev[k] != 0 and e in (1, 4)):
                            cl, cr = cr, cl
                        cur_edge_list[i] = np.array([n1, n2, cl, cr], int)

                    '''Flush non-adjacent edges into MSH file'''
                    next_edge_first = cur_edge_first + len(cur_edge_list)
                    msh.add_section(XF_Comment("Boundary edge {}:".format(be)))
                    zone_idx += 1
                    msh.add_section(XF_Face(zone_idx, cur_edge_first, next_edge_first - 1, bc_list[k][e - 1], FaceType.Linear, cur_edge_list))
                    msh.add_blank()
                    cur_edge_first = next_edge_first

        return msh

    @staticmethod
    def cell_idx_3d(coordinate, quadrant, dimension):
        i, j, k = coordinate
        u, v, w = dimension
        rdx = 1 + (u - 1) * (v - 1) * k + (u - 1) * j + i

        if quadrant == 1:
            return rdx
        elif quadrant == 2:
            return rdx - 1
        elif quadrant == 3:
            return rdx - u
        elif quadrant == 4:
            return rdx - (u - 1)
        elif quadrant == 5:
            return rdx - (u - 1) * (v - 1)
        elif quadrant == 6:
            return rdx - (u - 1) * (v - 1) - 1
        elif quadrant == 7:
            return rdx - (u - 1) * (v - 1) - (u - 1) + 1
        elif quadrant == 8:
            return rdx - (u - 1) * (v - 1) - (u - 1)
        else:
            raise ValueError("Invalid quadrant index.")

    @staticmethod
    def pnt_idx_3d(coordinate, norm_dir, dimension):
        if norm_dir not in ('X', 'x', 'Y', 'y', 'Z', 'z', 'U', 'u', 'V', 'v', 'W', 'w'):
            raise ValueError("Invalid norm direction representation.")

        i, j, k = coordinate
        u, v, w = dimension
        t = (u + 1) * (v + 1) * k + (u + 1) * j + i
        if norm_dir in ('X', 'x', 'U', 'u'):
            return t, t + (v + 1), t + (v + 2), t + 1
        elif norm_dir in ('Y', 'y', 'V', 'v'):
            return t, t + (u + 1), t + (u + 2), t + 1
        else:
            return t, t + (u + 1), t + u, t - 1

    @classmethod
    def from_str3d(cls, grid, bc=(BCType.VelocityInlet, BCType.Outflow, BCType.Wall, BCType.Wall, BCType.Wall, BCType.Wall)):
        """
        将3维单块结构网格转化为Fluent MSH格式
        :param grid: 结构网格
        :param bc: 6个面的边界条件
        :return: 可用于后续生成msh文件的XF_MSH对象
        :rtype: XF_MSH
        """

        '''Basic variables'''
        u, v, w = grid.shape[:3]
        dimension = 3
        zone_idx = 0
        total_cell = (u - 1) * (v - 1) * (w - 1)
        total_pnt = u * v * w
        total_face = u * (v - 1) * (w - 1) + v * (w - 1) * (u - 1) + w * (u - 1) * (v - 1)

        '''Initialize MSH file'''
        msh = cls()
        msh.add_section(XF_Header())
        msh.add_blank()
        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(dimension))
        msh.add_blank()

        '''Flush cell info to MSH file'''
        msh.add_section(XF_Comment("Cell Declaration:"))
        msh.add_section(XF_Cell.declaration(total_cell))
        msh.add_blank()
        zone_idx += 1
        msh.add_section(XF_Comment("Cell Info:"))
        msh.add_section(XF_Cell(zone_idx, 1, total_cell, CellType.Fluid, CellElement.Hexahedral))
        msh.add_blank()

        '''Collect point coordinates'''
        pnt_idx = 0
        pnt_list = np.empty((total_pnt, dimension), float)

        for k in range(w):
            for j in range(v):
                for i in range(u):
                    cls.dimensional_copy(pnt_list[pnt_idx], grid[i][j][k], dimension)

        '''Flush pnt info to MSH file'''
        msh.add_section(XF_Comment("Point Declaration:"))
        msh.add_section(XF_Node.declaration(total_pnt))
        msh.add_blank()
        zone_idx += 1
        msh.add_section(XF_Comment("Point Coordinates:"))
        msh.add_section(XF_Node(zone_idx, 1, total_pnt, NodeType.Any, dimension, pnt_list))
        msh.add_blank()

        '''Build interior adjacent description'''
        ifn = (u - 1) * (v - 1) * (w - 2) + (v - 1) * (w - 1) * (u - 2) + (w - 1) * (u - 1) * (v - 2)
        inter_face_desc = np.empty((ifn, 4), int)
        inter_face_cnt = 0

        '''Flush face declaration to MSH file'''
        msh.add_section(XF_Comment("Face Declaration:"))
        msh.add_section(XF_Face.declaration(total_face))
        msh.add_blank()

    @classmethod
    def from_str3d_multi(cls, grid_list, bc_list, adj_info):
        pass
