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

        if di != 0 and dj != 0:
            raise AssertionError("Invalid coordinates!")

        if di == 0:
            if dj == 1:
                return XF_MSH.cell_index(i, j, 2, U, V), XF_MSH.cell_index(i, j, 1, U, V)
            elif dj == -1:
                return XF_MSH.cell_index(i, j, 4, U, V), XF_MSH.cell_index(i, j, 3, U, V)
            else:
                raise ValueError("Invalid coordinates!")
        else:
            if di == 1:
                return XF_MSH.cell_index(i, j, 1, U, V), XF_MSH.cell_index(i, j, 4, U, V)
            elif di == -1:
                return XF_MSH.cell_index(i, j, 3, U, V), XF_MSH.cell_index(i, j, 2, U, V)
            else:
                raise ValueError("Invalid coordinates!")

    @staticmethod
    def dimensional_copy(dst, src, dim):
        for i in range(dim):
            dst[i] = src[i]

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
        :param adj_info: 依次每个单块结构网格的邻接关系序列
        :return: 可用于后续生成msh文件的XF_MSH对象
        :rtype: XF_MSH
        """

        '''Set up basic variables'''
        dimension = 2
        blk_num = len(blk_list)
        blk_shape = np.empty((blk_num, 2), int)
        adj_desc = np.zeros((blk_num, 4, 3), int)

        for k, blk in enumerate(blk_list):
            blk_shape[k] = blk.shape[:2]

        for entry in adj_info:
            b1, e1 = entry[0]
            b2, e2 = entry[1]
            reverse = 1 if entry[2] else 0

            adj_desc[b1][e1 - 1] = np.array([b2, e2, reverse], int)
            adj_desc[b2][e2 - 1] = np.array([b1, e1, reverse], int)

        '''Initialize MSH file'''
        msh = cls()
        zone_idx = 0
        msh.add_section(XF_Header())
        msh.add_blank()

        msh.add_section(XF_Comment("Dimension:"))
        msh.add_section(XF_Dimension(dimension))
        msh.add_blank()

        '''Counting cells'''
        cell_start = np.zeros(blk_num + 1, int)

        def cell_num(u, v):
            return (u - 1) * (v - 1)

        for k in range(blk_num):
            cur_cell_cnt = cell_num(blk.shape[0], blk.shape[1])
            cell_start[k + 1] = cell_start[k] + cur_cell_cnt

        '''Flush cell info to MSH file'''
        msh.add_section(XF_Comment("Cell part:"))
        msh.add_section(XF_Cell.declaration(cell_start[-1]))
        zone_idx += 1
        msh.add_section(XF_Cell(zone_idx, 1, cell_start[-1], CellType.Fluid, CellElement.Quadrilateral))
        msh.add_blank()

        '''Counting points'''
        pnt_num = 0
        for blk in blk_list:
            pnt_num += blk.shape[0] * blk.shape[1]

        def blk_edge_pnt(blk_idx, edge):
            return blk_shape[blk_idx][0] if edge in (1, 3) else blk_shape[blk_idx][1]

        for entry in adj_info:
            b1, e1 = entry[0]
            cur_pnt_num = blk_edge_pnt(b1, e1)
            pnt_num -= cur_pnt_num

        def boundary_pnt_num(u, v):
            return 2 * (u + v - 2)

        def boundary_coordinate(idx, u, v):
            s1, s2, s3, s4 = u - 1, u + v - 2, 2 * u + v - 3, boundary_pnt_num(u, v)
            if idx < s1:
                return idx, 0
            elif idx < s2:
                return s1, idx - s1
            elif idx < s3:
                return s3 - idx, v - 1
            elif idx < s4:
                return 0, s4 - idx
            else:
                raise ValueError("Boundary point {} goes beyond maximum index of current block".format(idx))

        def boundary_pnt_idx(i, j, u, v):
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

        bc_flag = []
        for blk in blk_list:
            cu = blk.shape[0]
            cv = blk.shape[1]
            cur_pnt_num = boundary_pnt_num(cu, cv)
            flag_ary = np.full(cur_pnt_num, -1, int)
            bc_flag.append(flag_ary)

        pnt_list = np.zeros((pnt_num, dimension))
        pnt_idx = 0

        def get_edge_idx(i, j, u, v):
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

        def get_counterpart_idx(k1, e1, i, j, k2, e2, r):
            u1, v1 = blk_shape[k1]
            u2, v2 = blk_shape[k2]

            if e1 in (1, 3):
                t = u1 - 1 - i if r else i
            elif e1 in (2, 4):
                t = v1 - 1 - j if r else j
            else:
                raise ValueError('Invalid edge index: {}'.format(e1))

            ii = 0
            jj = 0
            if e2 == 1:
                ii = t
            elif e2 == 2:
                jj = t
            elif e2 == 3:
                ii = t
                jj = v2 - 1
            elif e2 == 4:
                ii = u2 - 1
                jj = t
            else:
                raise ValueError("Invalid counterpart edge index: {}".format(e2))

            return boundary_pnt_idx(ii, jj, u2, v2)

        '''Handling internal points first'''
        internal_pnt_start = np.zeros(blk_num + 1, int)
        for k, blk in enumerate(blk_list):
            cu = blk.shape[0]
            cv = blk.shape[1]
            for j in range(1, cv - 1):
                for i in range(1, cu - 1):
                    cls.dimensional_copy(pnt_list[pnt_idx], blk[i][j], dimension)
                    pnt_idx += 1
            internal_pnt_start[k + 1] = internal_pnt_start[k] + pnt_idx

        '''Handling boundary points afterwards'''
        for k, blk in enumerate(blk_list):
            cu = blk.shape[0]
            cv = blk.shape[1]
            cn = boundary_pnt_num(cu, cv)
            for w in range(cn):
                if bc_flag[k][w] != -1:  # Has been marked
                    continue

                i, j = boundary_coordinate(w, cu, cv)
                edge_idx = get_edge_idx(i, j, cu, cv)

                '''Check if it's an interior pnt'''
                is_interior = False
                for e in edge_idx:
                    if adj_desc[k][e - 1][1] != 0:
                        is_interior = True
                        break

                if not is_interior:
                    cls.dimensional_copy(pnt_list[pnt_idx], blk[i][j], dimension)
                    pnt_idx += 1
                else:
                    '''Inspect neighbours'''
                    has_assigned = False
                    cur_adj = []
                    for e in edge_idx:
                        k2, e2, r = adj_desc[k][e - 1]
                        cp_idx = get_counterpart_idx(k, e, i, j, k2, e2, r)
                        cur_adj.append((k2, cp_idx))
                        if bc_flag[k2][cp_idx] != -1:
                            bc_flag[k][w] = bc_flag[k2][cp_idx]
                            has_assigned = True
                            break

                    if not has_assigned:
                        cls.dimensional_copy(pnt_list[pnt_idx], blk[i][j], dimension)
                        bc_flag[k][w] = pnt_idx
                        for e in cur_adj:
                            bc_flag[e[0]][e[1]] = pnt_idx
                        pnt_idx += 1

        '''Flush cell info to MSH file'''
        msh.add_section(XF_Comment("Point part:"))
        msh.add_section(XF_Node.declaration(pnt_num, dimension))
        zone_idx += 1
        msh.add_section(XF_Node(zone_idx, 1, pnt_num, NodeType.Any, dimension, pnt_list))
        msh.add_blank()

        '''Counting edges'''
        edge_num = 0

        def blk_edge_num(u, v):
            return (u - 1) * v + (v - 1) * u

        def boundary_edge_num(u, v, e=None):
            if e is None:
                return 2 * (u + v) - 4
            elif e in (1, 3):
                return u - 1
            elif e in (2, 4):
                return v - 1
            else:
                raise ValueError("Invalid edge index when counting boundary edges, should be in (1, 2, 3, 4) but get {}".format(e))

        for k in range(blk_num):
            edge_num += blk_edge_num(blk_shape[k][0], blk_shape[1])

        for entry in adj_info:
            k1, e1 = entry[0]
            cu, cv = blk_shape[k1]
            edge_num -= boundary_edge_num(cu, cv, e1)

        '''Flush edge declaration to MSH file'''
        msh.add_section(XF_Comment("Edge part:"))
        msh.add_section(XF_Face.declaration(edge_num))

        '''Handle internal edges in each blk'''
        for k, blk in enumerate(blk_list):
            cu, cv = blk_shape[k]
            pass
