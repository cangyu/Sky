import numpy as np
from abc import abstractmethod
from enum import Enum, unique
import platform


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

    @abstractmethod
    def write(self, out_stream):
        pass


class XF_Comment(XF_Section):
    def __init__(self, msg=''):
        """
        注释
        :param msg: 注释信息 
        :type msg: str
        """

        super(XF_Comment, self).__init__(0)
        self.msg = msg

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.msg)

    def write(self, out_stream):
        out_stream.write("({} \"{}\")".format(self.index, self.msg))


class XF_Header(XF_Section):
    def __init__(self, info=''):
        """
        To identify the program that wrote the file.
        :param info: Header info.
        :type info: str
        """

        super(XF_Header, self).__init__(1)
        self.info = "Grid generated with Python " + platform.python_version() if info == '' else info

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.info)

    def write(self, out_stream):
        out_stream.write("({} \"{}\")".format(self.index, self.info))


class XF_Dimension(XF_Section):
    def __init__(self, dim):
        """
        Dimension info of the grid.
        :param dim: Dimension 
        :type dim: int
        """

        super(XF_Dimension, self).__init__(2)
        self.ND = dim

    def build_content(self):
        self.formatted_content = "({} {})".format(self.index, self.ND)

    def write(self, out_stream):
        out_stream.write("({} {})".format(self.index, self.ND))


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

    def write(self, out_stream):
        if self.pts is None:
            out_stream.write("({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.node_type)[2:]))
        else:
            out_stream.write("({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.node_type)[2:], self.ND))
            n, dim = self.pts.shape
            for i in range(n):
                out_stream.write("\n{}".format(self.pts[i][0]))
                for d in range(1, self.ND):
                    out_stream.write(" {}".format(self.pts[i][d]))
            out_stream.write("))")


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

    def write(self, out_stream):
        if self.face_info is not None:
            out_stream.write("({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.bc_type)[2:], hex(self.face_type)[2:]))
            for fc in self.face_info:
                out_stream.write("\n{}".format(hex(fc[0])[2:]))
                cfn = len(fc)
                for cfi in range(1, cfn):
                    out_stream.write(" {}".format(hex(fc[cfi])[2:]))
            out_stream.write('))')
        else:
            out_stream.write("({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.face_type)[2:]))


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

    def write(self, out_stream):
        if self.cell_info is not None:
            out_stream.write("({} ({} {} {} {} {})(".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:], hex(self.elem_type)[2:]))
            out_stream.write("\n{}".format(hex(self.cell_info[0])[2:]))
            cn = len(self.cell_info)
            for ci in range(1, cn):
                out_stream.write(" {}".format(hex(self.cell_info[ci])[2:]))
            out_stream.write("))")
        else:
            if self.cell_type == 0:
                out_stream.write("({} ({} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:]))  # declaration
            else:
                out_stream.write("({} ({} {} {} {} {}))".format(self.index, hex(self.zone_id)[2:], hex(self.first_index)[2:], hex(self.last_index)[2:], hex(self.cell_type)[2:], hex(self.elem_type)[2:]))


class XF_MSH(object):
    def __init__(self):
        """
        ANSYS Fluent MSH文件
        """

        self.section_list = []

    def add_section(self, section):
        """
        向MSH文件添加信息段
        :param section: 信息段
        :type section: XF_Section
        :return: None
        """

        self.section_list.append(section)

    def save(self, fn):
        """
        输出MSH文件
        :param fn: 输出文件名
        :type fn: str
        :return: None
        """

        print("Writing ANSYS Fluent MSH file: \'{}\' ...".format(fn))
        msh = open(fn, 'w')
        for sec in self.section_list:
            sec.write(msh)
            msh.write("\n")
        msh.close()
        print("MSH file \'{}\' output done!".format(fn))

    @staticmethod
    def pnt_idx_2d(i, j, u):
        """
        坐标对应的序号(Starting from 1)
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
    def pnt_idx_3d(coordinate, dimension):
        """
        3维单块Block中某点的序号
        :param coordinate: 点的3维坐标
        :param dimension: Block的尺寸
        :return: 该点按照(x > y > z)的优先级排序下的序号，从1开始
        :rtype: int
        """

        i, j, k = coordinate
        u, v, w = dimension

        if i >= u or j >= v or k >= w:
            raise AssertionError("Invalid coordinate, boundary exceeded.")

        return 1 + u * v * k + u * j + i

    @staticmethod
    def pnt_idx_list_3d(coordinate, norm_dir, dimension):
        """
        求一点周围4个点的序号
        :param coordinate: 指定点的3维坐标(i, j, k)
        :param norm_dir: Face的法方向
        :param dimension: 3维Block的尺寸(u, v, w)
        :return: 从起始点coordinate开始，拇指指向norm_dir, 按右手法则绕一圈依次经过的4个点的序号(Starting from 1)
        """

        if norm_dir not in ('X', 'x', 'U', 'u', 'Y', 'y', 'V', 'v', 'Z', 'z', 'W', 'w'):
            raise AssertionError("Invalid norm direction representation.")

        i, j, k = coordinate
        if norm_dir in ('X', 'x', 'U', 'u'):
            t1 = XF_MSH.pnt_idx_3d(coordinate, dimension)
            t2 = XF_MSH.pnt_idx_3d((i, j + 1, k), dimension)
            t3 = XF_MSH.pnt_idx_3d((i, j + 1, k + 1), dimension)
            t4 = XF_MSH.pnt_idx_3d((i, j, k + 1), dimension)
            return t1, t2, t3, t4
        elif norm_dir in ('Y', 'y', 'V', 'v'):
            t1 = XF_MSH.pnt_idx_3d(coordinate, dimension)
            t2 = XF_MSH.pnt_idx_3d((i, j, k + 1), dimension)
            t3 = XF_MSH.pnt_idx_3d((i + 1, j, k + 1), dimension)
            t4 = XF_MSH.pnt_idx_3d((i + 1, j, k), dimension)
            return t1, t2, t3, t4
        else:
            t1 = XF_MSH.pnt_idx_3d(coordinate, dimension)
            t2 = XF_MSH.pnt_idx_3d((i + 1, j, k), dimension)
            t3 = XF_MSH.pnt_idx_3d((i + 1, j + 1, k), dimension)
            t4 = XF_MSH.pnt_idx_3d((i, j + 1, k), dimension)
            return t1, t2, t3, t4

    @staticmethod
    def cell_idx_quadrant_2d(i, j, d, u, v):
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
    def cell_idx_3d(coordinate, dimension):
        """
        计算3维Block中指定Cell的序号
        :param coordinate: 原点坐标
        :param dimension: Block尺寸
        :return: 以给定点为原点，正方向沿用Block的正方向，位于第1象限处的Cell的序号
        :rtype: int
        """

        i, j, k = coordinate
        u, v, w = dimension

        if i >= u - 1 or j >= v - 1 or k >= w - 1:
            raise AssertionError("Cell not Exist, boundary exceeded.")

        return 1 + (u - 1) * (v - 1) * k + (u - 1) * j + i

    @staticmethod
    def cell_idx_quadrant_3d(coordinate, quadrant, dimension):
        """
        计算3维Block中指定Cell的序号
        :param coordinate: 原点坐标
        :param quadrant: 象限号
        :type quadrant: int
        :param dimension: Block尺寸
        :return: 以给定点为原点，正方向沿用Block的正方向，位于指定象限处的Cell的序号
        :rtype: int
        """

        if quadrant < 1 or quadrant > 8:
            raise AssertionError("Invalid quadrant index.")

        i, j, k = coordinate

        if quadrant == 1:
            return XF_MSH.cell_idx_3d(coordinate, dimension)
        elif quadrant == 2:
            return XF_MSH.cell_idx_3d((i - 1, j, k), dimension)
        elif quadrant == 3:
            return XF_MSH.cell_idx_3d((i - 1, j - 1, k), dimension)
        elif quadrant == 4:
            return XF_MSH.cell_idx_3d((i, j - 1, k), dimension)
        elif quadrant == 5:
            return XF_MSH.cell_idx_3d((i, j, k - 1), dimension)
        elif quadrant == 6:
            return XF_MSH.cell_idx_3d((i - 1, j, k - 1), dimension)
        elif quadrant == 7:
            return XF_MSH.cell_idx_3d((i - 1, j - 1, k - 1), dimension)
        else:
            return XF_MSH.cell_idx_3d((i, j - 1, k - 1), dimension)

    @staticmethod
    def intersect_2d(p1, p2, u, v):
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
                return XF_MSH.cell_idx_quadrant_2d(i, j, 2, u, v), XF_MSH.cell_idx_quadrant_2d(i, j, 1, u, v)
            elif dj == -1:
                return XF_MSH.cell_idx_quadrant_2d(i, j, 4, u, v), XF_MSH.cell_idx_quadrant_2d(i, j, 3, u, v)
            else:
                raise ValueError("Invalid coordinates!")
        else:
            if di == 1:
                return XF_MSH.cell_idx_quadrant_2d(i, j, 1, u, v), XF_MSH.cell_idx_quadrant_2d(i, j, 4, u, v)
            elif di == -1:
                return XF_MSH.cell_idx_quadrant_2d(i, j, 3, u, v), XF_MSH.cell_idx_quadrant_2d(i, j, 2, u, v)
            else:
                raise ValueError("Invalid coordinates!")

    @classmethod
    def from_str2d(cls, grid, bc=(BCType.Wall, BCType.VelocityInlet, BCType.Wall, BCType.Outflow)):
        """
        从2维结构网格构建Fluent msh文件
        :param grid: 2D structural grid
        :param bc: Boundary Type list
        :return: 可用于后续生成msh文件的XF_MSH对象
        :rtype: XF_MSH
        """

        '''Basic variables'''
        cu, cv = grid.shape[:2]
        dim = 2
        total_node = node_num(cu, cv)
        total_cell = cell_num(cu, cv)
        total_edge = face_num(cu, cv)

        '''Initialize MSH file'''
        msh = cls()
        msh.add_section(XF_Header())
        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(dim))
        msh.add_section(XF_Comment("Declaration:"))
        msh.add_section(XF_Cell.declaration(total_cell))
        msh.add_section(XF_Face.declaration(total_edge))
        msh.add_section(XF_Node.declaration(total_node))
        zone_idx = 0

        '''Flush Cell info into MSH file'''
        msh.add_section(XF_Comment("Cells:"))
        zone_idx += 1
        msh.add_section(XF_Cell(zone_idx, 1, total_cell, CellType.Fluid, CellElement.Quadrilateral))

        '''Assemble node list'''
        tk = 0
        node_list = np.empty((total_node, dim), float)
        for j in range(cv):
            for i in range(cu):
                dimensional_copy(node_list[tk], grid[i][j], dim)  # 节点序号规定为沿x方向优先，依次递增
                tk += 1

        '''Flush Node info into MSH file'''
        msh.add_section(XF_Comment("Nodes:"))
        zone_idx += 1
        msh.add_section(XF_Node(zone_idx, 1, total_node, NodeType.Any, dim, node_list))

        '''Face'''
        fc = 0
        cur_face_num = cu - 1
        face_list = np.empty((cur_face_num, 4), int)
        for i in range(cur_face_num):
            face_list[i][0] = cls.pnt_idx_2d(i, 0, cu)
            face_list[i][1] = cls.pnt_idx_2d(i + 1, 0, cu)
            cl, cr = cls.intersect_2d((i, 0), (i + 1, 0), cu, cv)
            face_list[i][2] = cl
            face_list[i][3] = cr
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + cur_face_num, bc[0], FaceType.Linear, face_list))  # C1
        fc += cur_face_num

        cur_face_num = cv - 1
        face_list = np.empty((cur_face_num, 4), int)
        for j in range(cur_face_num):
            face_list[j][0] = cls.pnt_idx_2d(0, j, cu)
            face_list[j][1] = cls.pnt_idx_2d(0, j + 1, cu)
            cl, cr = cls.intersect_2d((0, j), (0, j + 1), cu, cv)
            face_list[j][2] = cl
            face_list[j][3] = cr
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + cur_face_num, bc[1], FaceType.Linear, face_list))  # C2
        fc += cur_face_num

        cur_face_num = cu - 1
        face_list = np.empty((cur_face_num, 4), int)
        for i in range(cur_face_num):
            face_list[i][0] = cls.pnt_idx_2d(i, cv - 1, cu)
            face_list[i][1] = cls.pnt_idx_2d(i + 1, cv - 1, cu)
            cl, cr = cls.intersect_2d((i, cv - 1), (i + 1, cv - 1), cu, cv)
            face_list[i][2] = cl
            face_list[i][3] = cr
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + cur_face_num, bc[2], FaceType.Linear, face_list))  # C3
        fc += cur_face_num

        cur_face_num = cv - 1
        face_list = np.empty((cur_face_num, 4), int)
        for j in range(cur_face_num):
            face_list[j][0] = cls.pnt_idx_2d(cu - 1, j, cu)
            face_list[j][1] = cls.pnt_idx_2d(cu - 1, j + 1, cu)
            cl, cr = cls.intersect_2d((cu - 1, j), (cu - 1, j + 1), cu, cv)
            face_list[j][2] = cl
            face_list[j][3] = cr
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + cur_face_num, bc[3], FaceType.Linear, face_list))  # C4
        fc += cur_face_num

        cur_face_num = total_edge - fc
        face_list = np.empty((cur_face_num, 4), int)
        tk = 0
        for i in range(1, cu - 1):
            for j in range(cv - 1):
                face_list[tk][0] = cls.pnt_idx_2d(i, j, cu)
                face_list[tk][1] = cls.pnt_idx_2d(i, j + 1, cu)
                cl, cr = cls.intersect_2d((i, j), (i, j + 1), cu, cv)
                face_list[tk][2] = cl
                face_list[tk][3] = cr
                tk += 1
        for j in range(1, cv - 1):
            for i in range(cu - 1):
                face_list[tk][0] = cls.pnt_idx_2d(i, j, cu)
                face_list[tk][1] = cls.pnt_idx_2d(i + 1, j, cu)
                cl, cr = cls.intersect_2d((i, j), (i + 1, j), cu, cv)
                face_list[tk][2] = cl
                face_list[tk][3] = cr
                tk += 1
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + cur_face_num, BCType.Interior, FaceType.Linear, face_list))  # Interior

        return msh

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
        dim = 3
        blk_shape = (u, v, w)
        total_cell = cell_num(u, v, w)
        total_pnt = node_num(u, v, w)
        total_face = face_num(u, v, w)

        '''Initialize MSH file and flush cell info'''
        msh = cls()
        msh.add_section(XF_Header())
        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(dim))
        msh.add_section(XF_Comment("Declaration:"))
        msh.add_section(XF_Cell.declaration(total_cell))
        msh.add_section(XF_Node.declaration(total_pnt))
        msh.add_section(XF_Face.declaration(total_face))
        zone_idx = 1
        msh.add_section(XF_Comment("Cell Info:"))
        msh.add_section(XF_Cell(zone_idx, 1, total_cell, CellType.Fluid, CellElement.Hexahedral))

        '''Collect point coordinates'''
        pnt_list = np.empty((total_pnt, dim), float)
        tk = 0
        for k in range(w):
            for j in range(v):
                for i in range(u):
                    dimensional_copy(pnt_list[tk], grid[i][j][k], dim)
                    tk += 1

        '''Flush pnt info to MSH file'''
        zone_idx += 1
        msh.add_section(XF_Comment("Point Coordinates:"))
        msh.add_section(XF_Node(zone_idx, 1, total_pnt, NodeType.Any, dim, pnt_list))

        '''Build interior adjacent description'''
        cur_face_num = internal_face_num(u, v, w)
        face_list = np.empty((cur_face_num, 6), int)
        fc = 0

        tk = 0
        for i in range(1, u - 1):
            for j in range(v - 1):
                for k in range(w - 1):
                    coord = (i, j, k)
                    t1, t2, t3, t4 = cls.pnt_idx_list_3d(coord, 'X', blk_shape)
                    face_list[tk][0] = t1
                    face_list[tk][1] = t2
                    face_list[tk][2] = t3
                    face_list[tk][3] = t4
                    face_list[tk][4] = XF_MSH.cell_idx_quadrant_3d(coord, 1, blk_shape)  # c0
                    face_list[tk][5] = XF_MSH.cell_idx_quadrant_3d(coord, 2, blk_shape)  # c1
                    tk += 1

        for j in range(1, v - 1):
            for k in range(w - 1):
                for i in range(u - 1):
                    coord = (i, j, k)
                    t1, t2, t3, t4 = cls.pnt_idx_list_3d(coord, 'Y', blk_shape)
                    face_list[tk][0] = t1
                    face_list[tk][1] = t2
                    face_list[tk][2] = t3
                    face_list[tk][3] = t4
                    face_list[tk][4] = XF_MSH.cell_idx_quadrant_3d(coord, 1, blk_shape)  # c0
                    face_list[tk][5] = XF_MSH.cell_idx_quadrant_3d(coord, 4, blk_shape)  # c1
                    tk += 1

        for k in range(1, w - 1):
            for i in range(u - 1):
                for j in range(v - 1):
                    coord = (i, j, k)
                    t1, t2, t3, t4 = cls.pnt_idx_list_3d(coord, 'Z', blk_shape)
                    face_list[tk][0] = t1
                    face_list[tk][1] = t2
                    face_list[tk][2] = t3
                    face_list[tk][3] = t4
                    face_list[tk][4] = XF_MSH.cell_idx_quadrant_3d(coord, 1, blk_shape)  # c0
                    face_list[tk][5] = XF_MSH.cell_idx_quadrant_3d(coord, 5, blk_shape)  # c1
                    tk += 1

        '''Flush Interior face info into MSH file'''
        zone_idx += 1
        msh.add_section(XF_Comment("Interior faces:"))
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + tk, BCType.Interior, FaceType.Quadrilateral, face_list))
        fc += tk

        '''X-Boundary faces'''
        cur_face_num = (v - 1) * (w - 1)
        tk = 0
        x0_face_desc = np.empty((cur_face_num, 6), int)
        x1_face_desc = np.empty((cur_face_num, 6), int)
        for j in range(v - 1):
            for k in range(w - 1):
                c1 = (0, j, k)
                c2 = (u - 1, j, k)
                f1, f2, f3, f4 = cls.pnt_idx_list_3d(c1, 'X', blk_shape)
                t1, t2, t3, t4 = cls.pnt_idx_list_3d(c2, 'X', blk_shape)

                x0_face_desc[tk][0] = f1
                x0_face_desc[tk][1] = f2
                x0_face_desc[tk][2] = f3
                x0_face_desc[tk][3] = f4
                x0_face_desc[tk][4] = XF_MSH.cell_idx_quadrant_3d(c1, 1, blk_shape)  # c0
                x0_face_desc[tk][5] = 0  # c1

                x1_face_desc[tk][0] = t1
                x1_face_desc[tk][1] = t2
                x1_face_desc[tk][2] = t3
                x1_face_desc[tk][3] = t4
                x1_face_desc[tk][4] = 0  # c0
                x1_face_desc[tk][5] = XF_MSH.cell_idx_quadrant_3d(c2, 2, blk_shape)  # c1
                tk += 1

        msh.add_section(XF_Comment("X-Boundary faces:"))
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + tk, bc[0], FaceType.Quadrilateral, x0_face_desc))
        fc += tk
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + tk, bc[1], FaceType.Quadrilateral, x1_face_desc))
        fc += tk

        '''Y-Boundary faces'''
        cur_face_num = (w - 1) * (u - 1)
        tk = 0
        y0_face_desc = np.empty((cur_face_num, 6), int)
        y1_face_desc = np.empty((cur_face_num, 6), int)
        for k in range(w - 1):
            for i in range(u - 1):
                c1 = (i, 0, k)
                c2 = (i, v - 1, k)
                f1, f2, f3, f4 = cls.pnt_idx_list_3d(c1, 'Y', blk_shape)
                t1, t2, t3, t4 = cls.pnt_idx_list_3d(c2, 'Y', blk_shape)

                y0_face_desc[tk][0] = f1
                y0_face_desc[tk][1] = f2
                y0_face_desc[tk][2] = f3
                y0_face_desc[tk][3] = f4
                y0_face_desc[tk][4] = XF_MSH.cell_idx_quadrant_3d(c1, 1, blk_shape)  # c0
                y0_face_desc[tk][5] = 0  # c1

                y1_face_desc[tk][0] = t1
                y1_face_desc[tk][1] = t2
                y1_face_desc[tk][2] = t3
                y1_face_desc[tk][3] = t4
                y1_face_desc[tk][4] = 0  # c0
                y1_face_desc[tk][5] = XF_MSH.cell_idx_quadrant_3d(c2, 4, blk_shape)  # c1
                tk += 1

        msh.add_section(XF_Comment("Y-Boundary faces:"))
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + tk, bc[2], FaceType.Quadrilateral, y0_face_desc))
        fc += tk
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + tk, bc[3], FaceType.Quadrilateral, y1_face_desc))
        fc += tk

        '''Z-Boundary faces'''
        cur_face_num = (u - 1) * (v - 1)
        tk = 0
        z0_face_desc = np.empty((cur_face_num, 6), int)
        z1_face_desc = np.empty((cur_face_num, 6), int)
        for i in range(u - 1):
            for j in range(v - 1):
                c1 = (i, j, 0)
                c2 = (i, j, w - 1)
                f1, f2, f3, f4 = cls.pnt_idx_list_3d(c1, 'Z', blk_shape)
                t1, t2, t3, t4 = cls.pnt_idx_list_3d(c2, 'Z', blk_shape)

                z0_face_desc[tk][0] = f1
                z0_face_desc[tk][1] = f2
                z0_face_desc[tk][2] = f3
                z0_face_desc[tk][3] = f4
                z0_face_desc[tk][4] = XF_MSH.cell_idx_quadrant_3d(c1, 1, blk_shape)  # c0
                z0_face_desc[tk][5] = 0  # c1

                z1_face_desc[tk][0] = t1
                z1_face_desc[tk][1] = t2
                z1_face_desc[tk][2] = t3
                z1_face_desc[tk][3] = t4
                z1_face_desc[tk][4] = 0  # c0
                z1_face_desc[tk][5] = XF_MSH.cell_idx_quadrant_3d(c2, 5, blk_shape)  # c1
                tk += 1

        msh.add_section(XF_Comment("Z-Boundary faces:"))
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + tk, bc[4], FaceType.Quadrilateral, z0_face_desc))
        fc += tk
        zone_idx += 1
        msh.add_section(XF_Face(zone_idx, fc + 1, fc + tk, bc[5], FaceType.Quadrilateral, z1_face_desc))
        fc += tk

        return msh

    @classmethod
    def from_str2d_multi(cls, blk_list, bc_list, adj_info):
        """
        根据给定的邻接关系，将2维多块结构网格转换成非结构形式，并赋予给定的边界条件,
        类似于ICEM中的'Convert to unstructured mesh'功能
        :param blk_list: 2维单块结构网格序列
        :param bc_list: 依次包括每个单块结构网格的BC
        :param adj_info: 依次包括每个单块结构网格的邻接关系序列
        :return: 可用于后续生成msh文件的XF_MSH对象
        :rtype: XF_MSH
        """

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
            if is_boundary_pnt(blk_shape[k], (i, j)):
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

            s1, s2, s3, s4 = u - 1, u + v - 2, 2 * u + v - 3, boundary_node_num((u, v))

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

            cn = boundary_node_num(blk_shape[k1], e1) - 1
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
        dim = 2
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
        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(dim))

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
        zone_idx += 1
        msh.add_section(XF_Comment("Cell Info:"))
        msh.add_section(XF_Cell(zone_idx, 1, total_cell, CellType.Fluid, CellElement.Quadrilateral))

        '''Counting points'''
        total_pnt = 0
        for k in range(total_blk):
            cu, cv = blk_shape[k]
            total_pnt += node_num(cu, cv)

        for entry in adj_info:
            b1, e1 = entry[0]
            b2, e2 = entry[1]
            if e1 != 0 and e2 != 0:
                total_pnt -= boundary_node_num(blk_shape[b1], e1)

        '''Flush point declaration msg to MSH file'''
        msh.add_section(XF_Comment("Point Declaration:"))
        msh.add_section(XF_Node.declaration(total_pnt))

        '''Local storage for recording boundary point index'''
        bc_flag = []
        for k, blk in enumerate(blk_list):
            cur_pnt_num = boundary_node_num(blk_shape[k])
            flag_ary = np.full(cur_pnt_num, -1, int)  # -1表示该点未分配序号
            bc_flag.append(flag_ary)

        pnt_list = np.empty((total_pnt, dim), float)
        pnt_idx = 0

        '''Handling internal points first'''
        internal_pnt_start = np.empty(total_blk + 1, int)
        internal_pnt_start[0] = 1
        for k, blk in enumerate(blk_list):
            cu, cv = blk_shape[k]
            for j in range(1, cv - 1):
                for i in range(1, cu - 1):
                    dimensional_copy(pnt_list[pnt_idx], blk[i][j], dim)
                    pnt_idx += 1
            internal_pnt_start[k + 1] = pnt_idx + 1

        '''Handling boundary points afterwards'''
        for k, blk in enumerate(blk_list):
            cu, cv = blk_shape[k]
            cn = boundary_node_num(blk_shape[k])
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
                    dimensional_copy(pnt_list[pnt_idx], blk[i][j], dim)
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
                        dimensional_copy(pnt_list[pnt_idx], blk[i][j], dim)
                        pnt_idx += 1
                        bc_flag[k][w] = pnt_idx
                        for e in cur_adj:
                            bc_flag[e[0]][e[1]] = pnt_idx

        '''Flush point coordinates to MSH file'''
        zone_idx += 1
        msh.add_section(XF_Comment("Point Coordinates:"))
        msh.add_section(XF_Node(zone_idx, 1, total_pnt, NodeType.Any, dim, pnt_list))

        '''Counting edges'''
        total_edge = 0
        for k in range(total_blk):
            cu, cv = blk_shape[k]
            total_edge += face_num(cu, cv)

        for entry in adj_info:
            k1, e1 = entry[0]
            k2, e2 = entry[1]
            if e1 != 0 and e2 != 0:
                total_edge -= boundary_face_num(blk_shape[k1], e1)

        '''Flush edge declaration to MSH file'''
        msh.add_section(XF_Comment("Edge Declaration:"))
        msh.add_section(XF_Face.declaration(total_edge))

        '''Handle internal edges in each blk'''
        cur_edge_first = 1
        for k, blk in enumerate(blk_list):
            cu, cv = blk_shape[k]
            cn = internal_face_num(cu, cv)
            cur_inter_edge = np.empty((cn, 4), int)
            ce = 0
            ts = cell_start[k] - 1

            '''Horizontal'''
            for j in range(1, cv - 1):
                for i in range(cu - 1):
                    n1 = get_pnt_idx(i, j, k)
                    n2 = get_pnt_idx(i + 1, j, k)
                    cl, cr = XF_MSH.intersect_2d((i, j), (i + 1, j), cu, cv)
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
                    cl, cr = XF_MSH.intersect_2d((i, j), (i, j + 1), cu, cv)
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
            cur_edge_first = next_edge_first

        '''Handle non-adjacent boundary'''
        be = 0
        for k in range(total_blk):
            for e in range(1, 5):
                if adj_desc[k][e - 1][1] == 0:
                    be += 1
                    cur_edge_cell = boundary_cell_list(k, e)
                    cur_edge_pnt = boundary_pnt_idx_list(k, e)
                    cn = boundary_node_num(blk_shape[k], e) - 1
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
                    cur_edge_first = next_edge_first

        return msh

    @classmethod
    def from_str3d_multi(cls, blk_list, bc_list, adj_info):
        """
        根据给定的邻接关系，将3维多块结构网格转换成非结构形式，并赋予给定的边界条件,
        类似于ICEM中的'Convert to unstructured mesh'功能
        Note that the face index and direction are defined as follows:
        --------------------------------------------------------------
            face    const   primary coordinate  secondary coordinate
        --------------------------------------------------------------
             1       U_min          V                    W
             2       U_max          V                    W
             3       V_min          W                    U
             4       V_max          W                    U
             5       W_min          U                    U
             6       W_max          U                    U
        --------------------------------------------------------------
        This notation is different from NASA's neural map file definition, but similar.
        :param blk_list: 3维单块结构网格序列
        :param bc_list: 依次包括每个单块结构网格的BC
        :param adj_info: 依次包括每个单块结构网格的邻接关系序列
        :return: 可用于后续生成msh文件的XF_MSH对象
        :rtype: XF_MSH
        """

        def shell_pnt_caste(_nu, _nv, _nw):
            _ret = np.empty(6, int)
            _ret[0] = _ret[1] = _nw * _nv
            _ret[2] = _ret[3] = (_nu - 2) * _nw
            _ret[4] = _ret[5] = (_nu - 2) * (_nv - 2)
            for i in range(1, 6):
                _ret[i] += _ret[i - 1]

            return _ret

        def shell_pnt_idx_from_coord(_k, _crd):
            _nu, _nv, _nw = blk_shape[_k]
            _ci, _cj, _ck = _crd

            if _ci == 0:
                return _cj + _ck * _nv
            elif _ci == _nu - 1:
                return blk_shell_caste[_k][0] + _cj + _ck * _nv
            else:
                if _cj == 0:
                    return blk_shell_caste[_k][1] + _ck + (_ci - 1) * _nw
                elif _cj == _nv - 1:
                    return blk_shell_caste[_k][2] + _ck + (_ci - 1) * _nw
                else:
                    if _ck == 0:
                        return blk_shell_caste[_k][3] + (_ci - 1) + (_cj - 1) * (_nu - 2)
                    elif _ck == _nw - 1:
                        return blk_shell_caste[_k][4] + (_ci - 1) + (_cj - 1) * (_nu - 2)
                    else:
                        raise ValueError("Not a boundary pnt.")

        def shell_pnt_coord_from_idx(_k, _p):
            _nu, _nv, _nw = blk_shape[_k]
            if _p < blk_shell_caste[_k][0]:
                _ci = 0
                _cj = _p % _nv
                _ck = _p // _nv
            elif _p < blk_shell_caste[_k][1]:
                cp = _p - blk_shell_caste[_k][0]
                _ci = _nu - 1
                _cj = cp % _nv
                _ck = cp // _nv
            elif _p < blk_shell_caste[_k][2]:
                cp = _p - blk_shell_caste[_k][1]
                _ci = 1 + cp // _nw
                _cj = 0
                _ck = cp % _nw
            elif _p < blk_shell_caste[_k][3]:
                cp = _p - blk_shell_caste[_k][2]
                _ci = 1 + cp // _nw
                _cj = _nv - 1
                _ck = cp % _nw
            elif _p < blk_shell_caste[_k][4]:
                cp = _p - blk_shell_caste[_k][3]
                _ci = cp % (_nu - 2) + 1
                _cj = cp // (_nu - 2) + 1
                _ck = 0
            elif _p < blk_shell_caste[_k][5]:
                cp = _p - blk_shell_caste[_k][4]
                _ci = cp % (_nu - 2) + 1
                _cj = cp // (_nu - 2) + 1
                _ck = _nw - 1
            else:
                raise ValueError("Invalid pnt index.")

            return _ci, _cj, _ck

        def shell_pnt_face_idx(_k, coord):
            _nu, _nv, _nw = blk_shape[_k]
            _ci, _cj, _ck = coord

            ret = []
            if _ci == 0:
                ret.append(1)
            if _ci == _nu - 1:
                ret.append(2)
            if _cj == 0:
                ret.append(3)
            if _cj == _nv - 1:
                ret.append(4)
            if _ck == 0:
                ret.append(5)
            if _ck == _nw - 1:
                ret.append(6)

            return ret

        def dominant_idx(_f):
            if _f in (1, 2):
                return 1, 2, 0
            elif _f in (3, 4):
                return 2, 0, 1
            elif _f in (5, 6):
                return 0, 1, 2
            else:
                raise ValueError("Invalid face index.")

        def invariant_coord(_nu, _nv, _nw, _f):
            if _f in (1, 3, 5):
                return 0
            elif _f == 2:
                return _nu - 1
            elif _f == 4:
                return _nv - 1
            elif _f == 6:
                return _nw - 1
            else:
                raise ValueError("Invalid face index \'{}\'.".format(_f))

        def get_counterpart_pnt_coord(_f1, _b2, _f2, _crd, _swp):
            _nu2, _nv2, _nw2 = blk_shape[_b2]
            di1 = dominant_idx(_f1)
            di2 = dominant_idx(_f2)
            ret = np.empty(3, int)
            ret[di2[0]] = _crd[di1[0]]
            ret[di2[1]] = _crd[di1[1]]
            ret[di2[2]] = invariant_coord(_nu2, _nv2, _nw2, _f2)
            if _swp:
                ret[di2[0]], ret[di2[1]] = ret[di2[1]], ret[di2[0]]

            return ret

        def pnt_around(_crd, _norm_dir):
            _ci, _cj, _ck = _crd
            if _norm_dir in ('X', 'x', 'U', 'u'):
                return np.array([[_ci, _cj, _ck], [_ci, _cj + 1, _ck], [_ci, _cj + 1, _ck + 1], [_ci, _cj, _ck + 1]], int)
            elif _norm_dir in ('Y', 'y', 'V', 'v'):
                return np.array([[_ci, _cj, _ck], [_ci, _cj, _ck + 1], [_ci + 1, _cj, _ck + 1], [_ci + 1, _cj, _ck]], int)
            elif _norm_dir in ('Z', 'z', 'W', 'w'):
                return np.array([[_ci, _cj, _ck], [_ci + 1, _cj, _ck], [_ci + 1, _cj + 1, _ck], [_ci, _cj + 1, _ck]], int)
            else:
                raise ValueError("Invalid normal direction.")

        def pnt_idx(_k, _crd):
            _nu, _nv, _nw = blk_shape[_k]
            _ci, _cj, _ck = _crd

            if is_boundary_pnt(blk_shape[_k], _crd):
                _t = shell_pnt_idx_from_coord(_k, _crd)
                return shell_pnt_idx[_k][_t]
            else:
                base = inter_pnt_start[_k]
                off = XF_MSH.pnt_idx_3d((_ci - 1, _cj - 1, _ck - 1), (_nu - 2, _nv - 2, _nw - 2)) - 1
                return base + off

        def cell_idx(_k, _crd, _quadrant):
            base = cell_start[_k]
            off = XF_MSH.cell_idx_quadrant_3d(_crd, _quadrant, blk_shape[_k]) - 1
            return base + off

        '''Basic variables'''
        dim = 3
        total_cell = 0
        total_pnt = 0
        total_face = 0
        blk_num = len(blk_list)
        blk_shape = np.empty((blk_num, 3), int)
        blk_shell_caste = np.empty((blk_num, 6), int)
        adj_desc = np.zeros((blk_num, 6, 3), int)
        cell_start = np.empty(blk_num, int)

        print("Converting multi-block grid with {} block(s) into ANSYS Fluent MSH format ...".format(blk_num))

        for k, blk in enumerate(blk_list):
            blk_shape[k] = blk.shape[:3]
            u, v, w = blk_shape[k]
            cell_start[k] = total_cell + 1
            total_cell += cell_num(u, v, w)
            total_pnt += node_num(u, v, w)
            total_face += face_num(u, v, w)
            blk_shell_caste[k] = shell_pnt_caste(u, v, w)

        for entry in adj_info:
            b1, f1 = entry[0]
            b2, f2 = entry[1]
            if f1 == 0 and f2 == 0:
                raise ValueError("Invalid adjacent settings.")

            if f1 != 0 and f2 != 0:
                swap = 1 if entry[3] else 0
                adj_desc[b1][f1 - 1][0] = b2
                adj_desc[b1][f1 - 1][1] = f2
                adj_desc[b1][f1 - 1][2] = swap
                adj_desc[b2][f2 - 1][0] = b1
                adj_desc[b2][f2 - 1][1] = f1
                adj_desc[b2][f2 - 1][2] = swap
                total_pnt -= boundary_node_num(blk_shape[b1], f1)
                total_face -= boundary_face_num(blk_shape[b1], f1)

        def edge_list_on_face(_f):
            """
            Get the 4 edge index on specified face.
            :param _f: face index
            :type _f: int
            :return: Edge index list, starting from 1
            """

            if _f == 1:
                return 2, 9, 6, 12
            elif _f == 2:
                return 4, 10, 8, 11
            elif _f == 3:
                return 9, 1, 10, 5
            elif _f == 4:
                return 12, 3, 11, 7
            elif _f == 5:
                return 1, 2, 3, 4
            elif _f == 6:
                return 5, 6, 7, 8
            else:
                raise ValueError("Invalid face index: \'{}\'.".format(_f))

        '''node cnt fix'''
        sub_cnt = np.zeros((blk_num, 12), int)
        for k in range(1, blk_num):
            for entry in adj_info:
                b1, f1 = entry[0]
                b2, f2 = entry[1]
                if f1 != 0 and f2 != 0:
                    if b1 < k and b2 == k:
                        for e in edge_list_on_face(f2):
                            sub_cnt[b2][e - 1] += 1
                    if b2 < k and b1 == k:
                        for e in edge_list_on_face(f1):
                            sub_cnt[b1][e - 1] += 1

        def pnt_num_on_edge(_k, _e):
            """
            Get the number of points on specified edge on certain block.
            :param _k: Block index. (Starting from 0)
            :type _k: int
            :param _e: Edge index. (Starting from 1)
            :type _e: int
            :return: Number of points on that edge.
            :rtype: int
            """

            _nu, _nv, _nw = blk_shape[_k]
            if _e in (1, 5, 7, 3):
                return _nu
            elif _e in (2, 6, 8, 4):
                return _nv
            elif _e in (9, 10, 11, 12):
                return _nw
            else:
                raise ValueError("Invalid edge index: \'{}\'.".format(_e))

        for k in range(blk_num):
            for e in range(12):
                if sub_cnt[k][e] == 2:
                    total_pnt += pnt_num_on_edge(k, e + 1)

        print("MSH grid general info: {} nodes, {} faces, {} cells.".format(total_pnt, total_face, total_cell))

        '''Initialize MSH file'''
        msh = cls()
        msh.add_section(XF_Header())
        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(dim))

        print("Building cells ...")

        '''Flush cell info to MSH file'''
        msh.add_section(XF_Comment("Cell Declaration:"))
        msh.add_section(XF_Cell.declaration(total_cell))
        msh.add_section(XF_Comment("Cell Info:"))
        zone_idx = 1
        msh.add_section(XF_Cell(zone_idx, 1, total_cell, CellType.Fluid, CellElement.Hexahedral))

        '''Points'''
        print("Building nodes(2 steps) ...")
        pnt_desc = np.empty((total_pnt, 3), float)
        inter_pnt_start = np.empty(blk_num, int)
        pnt_cnt = 0

        print("(1/2) Building interior nodes ...")

        '''Interior points in each block'''
        for k, blk in enumerate(blk_list):
            u, v, w = blk_shape[k]
            inter_pnt_start[k] = pnt_cnt + 1
            for ck in range(1, w - 1):
                for cj in range(1, v - 1):
                    for ci in range(1, u - 1):
                        dimensional_copy(pnt_desc[pnt_cnt], blk[ci][cj][ck], dim)
                        pnt_cnt += 1

        print("(2/2) Building boundary nodes ...")

        shell_pnt_idx = []
        shell_pnt_num = np.empty(blk_num, int)
        for k in range(blk_num):
            pn = boundary_node_num(blk_shape[k])
            shell_pnt_num[k] = pn
            shell_pnt_idx.append(np.zeros(pn, int))

        for k, blk in enumerate(blk_list):
            pn = shell_pnt_num[k]
            for pc in range(pn):
                if shell_pnt_idx[k][pc] != 0:
                    continue

                '''Record new pnt coordinate without duplication'''
                ci, cj, ck = shell_pnt_coord_from_idx(k, pc)
                dimensional_copy(pnt_desc[pnt_cnt], blk[ci][cj][ck], dim)
                pnt_cnt += 1

                '''Find all adjacent point with BFS strategy'''
                adj_blk_pnt = [(k, pc)]
                t = 0
                while t < len(adj_blk_pnt):
                    k1, p1 = adj_blk_pnt[t]
                    crd = shell_pnt_coord_from_idx(k1, p1)
                    face = shell_pnt_face_idx(k1, crd)
                    for f1 in face:
                        if adj_desc[k1][f1 - 1][1] != 0:
                            b2, f2, swap = adj_desc[k1][f1 - 1]
                            cp_crd = get_counterpart_pnt_coord(f1, b2, f2, crd, swap)
                            p2 = shell_pnt_idx_from_coord(b2, cp_crd)
                            ca = (b2, p2)
                            if ca not in adj_blk_pnt:
                                adj_blk_pnt.append(ca)
                    t += 1

                '''Mark all adjacent points'''
                for entry in adj_blk_pnt:
                    shell_pnt_idx[entry[0]][entry[1]] = pnt_cnt

        '''Flush pnt info to MSH file'''
        msh.add_section(XF_Comment("Point Declaration:"))
        msh.add_section(XF_Node.declaration(total_pnt))
        zone_idx += 1
        msh.add_section(XF_Comment("Point Coordinates:"))
        msh.add_section(XF_Node(zone_idx, 1, total_pnt, NodeType.Any, dim, pnt_desc))

        print("Building faces(2 steps) ...")

        '''Flush face declaration to MSH file'''
        msh.add_section(XF_Comment("Face Declaration:"))
        msh.add_section(XF_Face.declaration(total_face))

        print("(1/2) Building interior faces ...")

        '''Interior Faces'''
        face_start = 1
        for k in range(blk_num):
            u, v, w = blk_shape[k]
            ifn = internal_face_num(u, v, w)
            inter_face_desc = np.empty((ifn, 6), int)
            cur_face_cnt = 0

            for nu in range(1, u - 1):
                for nv in range(v - 1):
                    for nw in range(w - 1):
                        cur_crd = (nu, nv, nw)
                        crd_list = pnt_around(cur_crd, 'X')
                        inter_face_desc[cur_face_cnt][0] = pnt_idx(k, crd_list[0])
                        inter_face_desc[cur_face_cnt][1] = pnt_idx(k, crd_list[1])
                        inter_face_desc[cur_face_cnt][2] = pnt_idx(k, crd_list[2])
                        inter_face_desc[cur_face_cnt][3] = pnt_idx(k, crd_list[3])
                        inter_face_desc[cur_face_cnt][4] = cell_idx(k, cur_crd, 1)  # c0
                        inter_face_desc[cur_face_cnt][5] = cell_idx(k, cur_crd, 2)  # c1
                        cur_face_cnt += 1

            for nv in range(1, v - 1):
                for nw in range(w - 1):
                    for nu in range(u - 1):
                        cur_crd = (nu, nv, nw)
                        crd_list = pnt_around(cur_crd, 'Y')
                        inter_face_desc[cur_face_cnt][0] = pnt_idx(k, crd_list[0])
                        inter_face_desc[cur_face_cnt][1] = pnt_idx(k, crd_list[1])
                        inter_face_desc[cur_face_cnt][2] = pnt_idx(k, crd_list[2])
                        inter_face_desc[cur_face_cnt][3] = pnt_idx(k, crd_list[3])
                        inter_face_desc[cur_face_cnt][4] = cell_idx(k, cur_crd, 1)  # c0
                        inter_face_desc[cur_face_cnt][5] = cell_idx(k, cur_crd, 4)  # c1
                        cur_face_cnt += 1

            for nw in range(1, w - 1):
                for nu in range(u - 1):
                    for nv in range(v - 1):
                        cur_crd = (nu, nv, nw)
                        crd_list = pnt_around(cur_crd, 'Z')
                        inter_face_desc[cur_face_cnt][0] = pnt_idx(k, crd_list[0])
                        inter_face_desc[cur_face_cnt][1] = pnt_idx(k, crd_list[1])
                        inter_face_desc[cur_face_cnt][2] = pnt_idx(k, crd_list[2])
                        inter_face_desc[cur_face_cnt][3] = pnt_idx(k, crd_list[3])
                        inter_face_desc[cur_face_cnt][4] = cell_idx(k, cur_crd, 1)  # c0
                        inter_face_desc[cur_face_cnt][5] = cell_idx(k, cur_crd, 5)  # c1
                        cur_face_cnt += 1

            '''Flush Interior face info into MSH file'''
            zone_idx += 1
            msh.add_section(XF_Comment("Block {} interior faces:".format(k + 1)))
            msh.add_section(XF_Face(zone_idx, face_start, face_start + cur_face_cnt - 1, BCType.Interior, FaceType.Quadrilateral, inter_face_desc))
            face_start += cur_face_cnt

        def face_coordinate_indented(_dim, _f):
            _bfn = boundary_face_num(_dim, _f)
            _ret = np.empty((_bfn, 3), int)
            _nu, _nv, _nw = _dim
            _pc = 0

            if _f in (1, 2):
                d1c = 0 if _f == 1 else _nu - 1
                for _j in range(_nv - 1):
                    for _k in range(_nw - 1):
                        _ret[_pc][0] = d1c
                        _ret[_pc][1] = _j
                        _ret[_pc][2] = _k
                        _pc += 1

            if _f in (3, 4):
                d2c = 0 if _f == 3 else _nv - 1
                for _k in range(_nw - 1):
                    for _i in range(_nu - 1):
                        _ret[_pc][0] = _i
                        _ret[_pc][1] = d2c
                        _ret[_pc][2] = _k
                        _pc += 1

            if _f in (5, 6):
                d3c = 0 if _f == 5 else _nw - 1
                for _i in range(_nu - 1):
                    for _j in range(_nv - 1):
                        _ret[_pc][0] = _i
                        _ret[_pc][1] = _j
                        _ret[_pc][2] = d3c
                        _pc += 1

            return _ret

        def face_norm_dir(_f):
            """
            将face编号转换为对应的法方向
            :param _f: Index of a face.
            :type _f: int
            :return: Corresponding norm direction.
            """

            if _f in (1, 2):
                return 'X'
            elif _f in (3, 4):
                return 'Y'
            elif _f in (5, 6):
                return 'Z'
            else:
                raise ValueError("Invalid face index.")

        def cell_quadrant_on_face(_f):
            """
            计算边界面上Cell的象限
            :param _f: Index of a face.
            :type _f: int
            :return:与边界边上的Face相邻的Cell相对于约定的原点所在的象限
            :rtype: int
            """

            if _f in (1, 3, 5):
                return 1
            elif _f == 2:
                return 2
            elif _f == 4:
                return 4
            elif _f == 6:
                return 5
            else:
                raise ValueError("Invalid face index.")

        def calc_boundary_face_adj_info(_b, _f, _left):
            """
            计算Block上指定面的邻接信息
            :param _b: Index of block. (Starting from 0)
            :type _b: int
            :param _f: Index of face. (Starting from 1)
            :type _f: int
            :param _left: Indicate if the face is on the left side of the intersection.
            :type _left: bool
            :return: Adjacent info used to generate ANSYS Fluent MSH file.
            """

            '''Get all coordinates'''
            _crd = face_coordinate_indented(blk_shape[_b], _f)
            _nd = face_norm_dir(_f)
            _qud = cell_quadrant_on_face(_f)

            '''Extend to surrounding 4 point index'''
            _bfn = len(_crd)
            _ret = np.empty((_bfn, 6), int)
            _tk = 0
            while _tk < _bfn:
                crd_ard = pnt_around(_crd[_tk], _nd)
                _ret[_tk][0] = pnt_idx(_b, crd_ard[0])
                _ret[_tk][1] = pnt_idx(_b, crd_ard[1])
                _ret[_tk][2] = pnt_idx(_b, crd_ard[2])
                _ret[_tk][3] = pnt_idx(_b, crd_ard[3])
                _ret[_tk][4] = 0
                _ret[_tk][5] = cell_idx(_b, _crd[_tk], _qud)
                if not _left:
                    _ret[_tk][4], _ret[_tk][5] = _ret[_tk][5], _ret[_tk][4]
                _tk += 1

            return _ret

        def calc_interior_face_adj_info(_b1, _f1, _b2, _f2, _ref, _swp):
            """
            计算两个交接面上的邻接信息
            :param _b1: Block index of the left face
            :type _b1: int
            :param _f1: Face index of the left face within its block
            :type _f1: int
            :param _b2: Block index of the right face
            :type _b2: int
            :param _f2: Face index of the right face within its block
            :type _f2: int
            :param _ref: Indicate which face is selected as reference, 0 means the left, 1 means the right
            :type _ref: int
            :param _swp: Indicate if the two faces' primary coordinates are in the same order
            :type _swp: bool
            :return: Adjacent info of this interior face.
            """

            _bfn = boundary_face_num(blk_shape[_b1], _f1)
            if _bfn != boundary_face_num(blk_shape[_b2], _f2):
                raise AssertionError("Invalid adjacent info.")

            '''Get all coordinates and corresponding denotational info'''
            _nd1 = face_norm_dir(_f1)
            _nd2 = face_norm_dir(_f2)
            _qud1 = cell_quadrant_on_face(_f1)
            _qud2 = cell_quadrant_on_face(_f2)
            _crd1 = face_coordinate_indented(blk_shape[_b1], _f1)
            _crd2 = face_coordinate_indented(blk_shape[_b2], _f2)
            _crd = _crd1 if _ref == 0 else _crd2
            _nd = _nd1 if _ref == 0 else _nd2
            _b = _b1 if _ref == 0 else _b2

            '''Calculate surrounding 4 point index'''
            _ret = np.empty((_bfn, 6), int)
            _pc = 0
            while _pc < _bfn:
                crd_ard = pnt_around(_crd[_pc], _nd)
                for c in range(4):
                    _ret[_pc][c] = pnt_idx(_b, crd_ard[c])

                if _ref == 0:
                    _ret[_pc][4] = cell_idx(_b1, _crd1[_pc], _qud1)
                    _cp_crd = get_counterpart_pnt_coord(_f1, _b2, _f2, _crd1[_pc], _swp)
                    _ret[_pc][5] = cell_idx(_b2, _cp_crd, _qud2)
                else:
                    _ret[_pc][5] = cell_idx(_b2, _crd2[_pc], _qud2)
                    _cp_crd = get_counterpart_pnt_coord(_f2, _b1, _f1, _crd2[_pc], _swp)
                    _ret[_pc][4] = cell_idx(_b1, _cp_crd, _qud1)

                _ret[_pc][4], _ret[_pc][5] = _ret[_pc][5], _ret[_pc][4]
                _pc += 1

            return _ret

        print("(2/2) Building boundary faces ...")

        '''Boundary Faces'''
        for k, entry in enumerate(adj_info):
            b1, f1 = entry[0]
            b2, f2 = entry[1]
            ref_idx = entry[2]
            swp = entry[3]

            if f1 == 0 or f2 == 0:
                '''Non-Adjacent faces'''
                on_left = f2 == 0
                boundary_face_desc = calc_boundary_face_adj_info(b1, f1, True) if on_left else calc_boundary_face_adj_info(b2, f2, False)
                bc = bc_list[b1][f1 - 1] if on_left else bc_list[b2][f2 - 1]
            else:
                '''Adjacent faces'''
                boundary_face_desc = calc_interior_face_adj_info(b1, f1, b2, f2, ref_idx, swp)
                bc = BCType.Interior

            cur_face_cnt = len(boundary_face_desc)
            if f1 != 0:
                msh.add_section(XF_Comment("Blk {}, Boundary {} faces:".format(b1 + 1, f1)))
            else:
                msh.add_section(XF_Comment("Blk {}, Boundary {} faces:".format(b2 + 1, f2)))
            zone_idx += 1
            msh.add_section(XF_Face(zone_idx, face_start, face_start + cur_face_cnt - 1, bc, FaceType.Quadrilateral, boundary_face_desc))
            face_start += cur_face_cnt

        print("Conversion done!")

        return msh


def boundary_face_num(_dim, _f=None):
    if len(_dim) == 2:
        _u, _v = _dim
        if _f in (1, 2, 3, 4):
            return _u - 1 if _f in (1, 3) else _v - 1
        elif _f is None:
            return face_num(_u, _v) - internal_face_num(_u, _v)
        else:
            raise ValueError("Invalid input.")
    elif len(_dim) == 3:
        _u, _v, _w = _dim
        if _f in (1, 2, 3, 4, 5, 6):
            if _f in (1, 2):
                return cell_num(_v, _w)
            elif _f in (3, 4):
                return cell_num(_w, _u)
            elif _f in (5, 6):
                return cell_num(_u, _v)
        elif _f is None:
            return face_num(_u, _v, _w) - internal_face_num(_u, _v, _w)
    else:
        raise ValueError("Invalid input.")





def boundary_node_num(_dim, _f=None):
    if len(_dim) == 2:
        _u, _v = _dim
        if _f is None:
            return node_num(_u, _v) - node_num(_u - 2, _v - 2)
        elif _f in (1, 3):
            return _u
        elif _f in (2, 4):
            return _v
        else:
            raise ValueError('Invalid face index.')
    elif len(_dim) == 3:
        _u, _v, _w = _dim
        if _f is None:
            return node_num(_u, _v, _w) - node_num(_u - 2, _v - 2, _w - 2)
        elif _f in (1, 2):
            return node_num(_v, _w)
        elif _f in (3, 4):
            return node_num(_w, _u)
        elif _f in (5, 6):
            return node_num(_u, _v)
        else:
            raise ValueError('Invalid face index.')
    else:
        raise ValueError("Invalid input.")








def dimensional_copy(dst, src, dim):
    for i in range(dim):
        dst[i] = src[i]
