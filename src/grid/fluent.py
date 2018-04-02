#!/usr/bin/env python
# -*- coding: utf-8 -*-

import platform
from abc import abstractmethod
from enum import Enum, unique
import numpy as np

"""
Implementation of the ANSYS Fluent MSH File Format.

Note:
Refer from the appendix of ANSYS Fluent 15.0 User Manual. 
"""


class XFSection(object):
    def __init__(self, idx):
        """
        Abstract class of sections in the ANSYS Fluent MSH File Format.
        :param idx: Category of the section.
        :type idx: int
        """

        self.index = idx
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


class XFComment(XFSection):
    def __init__(self, msg=''):
        """
        Comment Section in the ANSYS Fluent MSH File Format.
        :param msg: The comment message.
        :type msg: str
        """

        super(XFComment, self).__init__(0)
        self.msg = msg

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.msg)

    def write(self, out_stream):
        out_stream.write("({} \"{}\")".format(self.index, self.msg))


class XFHeader(XFSection):
    def __init__(self, info=''):
        """
        Header Section in the ANSYS Fluent MSH File Format.
        It's used to identify the program that wrote the file.
        :param info: Header info.
        :type info: str
        """

        super(XFHeader, self).__init__(1)
        self.info = "Grid generated with Python " + platform.python_version() if info == '' else info

    def build_content(self):
        self.formatted_content = "({} \"{}\")".format(self.index, self.info)

    def write(self, out_stream):
        out_stream.write("({} \"{}\")".format(self.index, self.info))


class XFDimension(XFSection):
    def __init__(self, dim):
        """
        Dimension Section in the ANSYS Fluent MSH File Format.
        It's used to specified the dimension of the grid.
        :param dim: Dimension
        :type dim: int
        """

        assert dim in (2, 3)
        super(XFDimension, self).__init__(2)
        self.ND = dim

    def build_content(self):
        self.formatted_content = "({} {})".format(self.index, self.ND)

    def write(self, out_stream):
        out_stream.write("({} {})".format(self.index, self.ND))


@unique
class NodeType(Enum):
    Virtual, Any, Boundary = range(3)


class XFNode(XFSection):
    def __init__(self, zone, first, last, tp, dim, pts=None):
        """
        Node Section in the ANSYS Fluent MSH File Format.
        :param zone: 区域号
        :type zone: int
        :param first: 网格点起始序号(Starting from 1)
        :type first: int
        :param last:  网格点终止序号
        :type last: int
        :param tp: 网格点类型(TGrid usage)
        :type tp: NodeType
        :param dim: Indicate the dimensionality of the node data.
        :type dim: int
        :param pts: 网格点数组，(last-first+1)个元素
        """

        super(XFNode, self).__init__(10)
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


class XFFace(XFSection):
    def __init__(self, zone, first, last, bct, ft, face_info=None):
        """
        Face Section in the ANSYS Fluent MSH File Format.
        :param zone: 区域号
        :type zone: int
        :param first: 起始序号(Starting from 1)
        :type first: int
        :param last: 终止序号
        :type last: int
        :param bct: 边界属性
        :type bct: BCType
        :param ft: 形状类别
        :type ft: FaceType
        :param face_info: 边界信息，一维数组(last-first+1)个元素，每个元素中包含了邻接关系
        """

        super(XFFace, self).__init__(13)
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


class XFCell(XFSection):
    def __init__(self, zone: int, first: int, last: int, ct: CellType, ce: CellElement, cell_info=None):
        """
        Cell Section in the ANSYS Fluent MSH File Format.
        :param zone: 区域号
        :param first: 起始序号(Starting from 1)
        :param last: 终止序号
        :param ct: Cell type.
        :param ce: Cell element.
        :param cell_info: 单元类型信息，一维数组，(last-first+1)个元素
        """

        super(XFCell, self).__init__(12)
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


def pnt_circle_x(pnt):
    # [i, j, k], [i, j + 1, k], [i, j + 1, k + 1], [i, j, k + 1]
    ret = np.array([pnt] * 4)
    ret[1][1] += 1
    ret[2][1] += 1
    ret[2][2] += 1
    ret[3][2] += 1
    return ret


def pnt_circle_y(pnt):
    # [i, j, k], [i, j, k + 1], [i + 1, j, k + 1], [i + 1, j, k]
    ret = np.array([pnt] * 4)
    ret[1][2] += 1
    ret[2][0] += 1
    ret[2][2] += 1
    ret[3][0] += 1
    return ret


def pnt_circle_z(pnt):
    # [i, j, k], [i + 1, j, k], [i + 1, j + 1, k], [i, j + 1, k]
    ret = np.array([pnt] * 4)
    ret[1][0] += 1
    ret[2][0] += 1
    ret[2][1] += 1
    ret[3][1] += 1
    return ret


def pnt_circle(pnt, norm_dir):
    """
    Calculate the surrounding 4 points from given pnt
    with specified norm direction.
    :param pnt: Origin.
    :param norm_dir: Norm direction.
                     0-(X, x, I, i, U, u)
                     1-(Y, y, J, j, V, v)
                     2-(Z, z, K, k, W, w)
    :type norm_dir: int
    :return: Surrounding 4 points.
    """

    if norm_dir == 0:
        return pnt_circle_x(pnt)
    elif norm_dir == 1:
        return pnt_circle_y(pnt)
    elif norm_dir == 2:
        return pnt_circle_z(pnt)
    else:
        raise ValueError('Invalid norm direction.')


def pnt_circle_h(pnt):
    # [i, j], [i+1, j]
    ret = np.array([pnt] * 2)
    ret[1][0] += 1
    return ret


def pnt_circle_v(pnt):
    # [i, j], [i, j+1]
    ret = np.array([pnt] * 2)
    ret[1][1] += 1
    return ret


def xf_calc_boundary_info(entry, nmf):
    """
    Calculate the boundary adjacent info in Fluent MSH format with given description.
    :param entry: Adjacent description.
    :type entry: NMFEntry
    :param nmf: Neutral Map File holding global info.
    :type nmf: NeutralMapFile
    :return: Boundary adjacent info in MSH format.
    """

    dim = nmf.dim
    elem_num = 6 if dim == 3 else 4
    n = entry.face_num
    ret = np.empty((n, elem_num), int)

    n = 0
    if dim == 3:
        if entry.Type == 'ONE_TO_ONE':
            for x1 in range(entry.pri_node_num() - 1):
                for x2 in range(entry.sec_node_num() - 1):
                    p1 = nmf.calc_real_pnt(entry, (x1, x2), 1)
                    p2 = nmf.calc_real_pnt(entry, (x1, x2), 2)
                    norm_dir = NeutralMapFile.FACE_INVARIANT[entry.F1]
                    crd_list = pnt_circle(p1, norm_dir)
                    for t in range(4):
                        ret[n][t] = nmf.calc_pnt_idx(entry.B1, crd_list[t])
                    c1 = nmf.calc_cell_idx(entry.B1, p1, NeutralMapFile.CELL_QUADRANT_ON_FACE[entry.F1])
                    c2 = nmf.calc_cell_idx(entry.B2, p2, NeutralMapFile.CELL_QUADRANT_ON_FACE[entry.F2])
                    ret[n][4] = c1
                    ret[n][5] = c2
                    n += 1
        else:
            for x1 in range(entry.pri_node_num - 1):
                for x2 in range(entry.sec_node_num - 1):
                    p = nmf.calc_real_pnt(entry, (x1, x2), 1)
                    norm_dir = NeutralMapFile.FACE_INVARIANT[entry.F1]
                    crd_list = pnt_circle(p, norm_dir)
                    for t in range(4):
                        ret[n][t] = nmf.calc_pnt_idx(entry.B1, crd_list[t])
                    c1 = nmf.calc_cell_idx(entry.B1, p, NeutralMapFile.CELL_QUADRANT_ON_FACE[entry.F1])
                    c2 = 0
                    ret[n][4] = c1
                    ret[n][5] = c2
                    n += 1
    else:
        pass

    return ret


class FluentMSH(object):
    NMF2MSH_BC_DICT = {'ONE_TO_ONE': BCType.Interior,
                       'WALL': BCType.Wall,
                       'Symmetry-X': BCType.Symmetry,
                       'Symmetry-Y': BCType.Symmetry,
                       'Symmetry-Z': BCType.Symmetry,
                       'FAR': BCType.PressureFarField}

    def __init__(self):
        """
        ANSYS Fluent MSH File
        """

        self.xf_section = []

    def add(self, section):
        """
        向MSH文件添加信息段
        :param section: 信息段
        :type section: XFSection
        :return: None
        """

        self.xf_section.append(section)

    def clear(self):
        self.xf_section.clear()

    @property
    def size(self):
        return len(self.xf_section)

    def save(self, fn):
        """
        输出MSH文件
        :param fn: 输出文件名
        :type fn: str
        :return: None
        """

        f_out = open(fn, 'w')
        for i in range(self.size):
            self.xf_section[i].write(f_out)
            f_out.write('\n')
        f_out.close()

    @classmethod
    def from_nmf(cls, nmf):
        """
        Construct ANSYS Fluent MSH file from Neutral Map File.
        :param nmf: Neutral mapping description.
        :type nmf: NeutralMapFile
        :return: Structured grid in ANSYS Fluent MSH File Format.
        :rtype: FluentMSH
        """

        msh = cls()
        dim = nmf.dim
        zone_idx = 1
        face_idx = 1
        cell_num = nmf.cell_num
        pnt_num = nmf.pnt_num
        face_num = nmf.face_num

        '''Declaration'''
        msh.add(XFHeader())
        msh.add(XFComment('Dimension:'))
        msh.add(XFDimension(dim))
        msh.add(XFComment("Cell Declaration:"))
        msh.add(XFCell.declaration(cell_num))
        msh.add(XFComment("Point Declaration:"))
        msh.add(XFNode.declaration(pnt_num))
        msh.add(XFComment("Face Declaration:"))
        msh.add(XFFace.declaration(face_num))

        '''Cell'''
        msh.add(XFComment("Cell Info:"))
        cell_type = CellElement.Hexahedral if dim == 3 else CellElement.Quadrilateral
        msh.add(XFCell(zone_idx, 1, cell_num, CellType.Fluid, cell_type))
        zone_idx += 1

        '''Point'''
        msh.add(XFComment("Point Coordinates:"))
        msh.add(XFNode(zone_idx, 1, pnt_num, NodeType.Any, dim, nmf.all_pnt))
        zone_idx += 1

        '''Interior Face'''
        if dim == 3:
            for b in range(nmf.blk_num):
                cur_face_total = blk_internal_face_num(nmf.blk[b].shape)
                face = np.empty((cur_face_total, 6), int)
                u, v, w, _ = nmf.blk[b].shape
                n = 0

                for i in range(1, u - 1):
                    for j in range(v - 1):
                        for k in range(w - 1):
                            p = np.array([i, j, k])
                            crd_list = pnt_circle_x(p)
                            for t in range(4):
                                face[n][t] = nmf.calc_pnt_idx(k, crd_list[t])
                            face[n][4] = nmf.calc_cell_idx(k, p, 1)  # c0
                            face[n][5] = nmf.calc_cell_idx(k, p, 2)  # c1
                            n += 1

                for j in range(1, v - 1):
                    for k in range(w - 1):
                        for i in range(u - 1):
                            p = np.array([i, j, k])
                            crd_list = pnt_circle_y(p)
                            for t in range(4):
                                face[n][t] = nmf.calc_pnt_idx(k, crd_list[t])
                            face[n][4] = nmf.calc_cell_idx(k, p, 1)  # c0
                            face[n][5] = nmf.calc_cell_idx(k, p, 4)  # c1
                            n += 1

                for k in range(1, w - 1):
                    for i in range(u - 1):
                        for j in range(v - 1):
                            p = np.array([i, j, k])
                            crd_list = pnt_circle_z(p)
                            for t in range(4):
                                face[n][t] = nmf.calc_pnt_idx(k, crd_list[t])
                            face[n][4] = nmf.calc_cell_idx(k, p, 1)  # c0
                            face[n][5] = nmf.calc_cell_idx(k, p, 5)  # c1
                            n += 1

                msh.add(XFComment("Blk-{} interior faces:".format(b + 1)))
                msh.add(XFFace(zone_idx, face_idx, face_idx + n - 1, BCType.Interior, FaceType.Quadrilateral, face))
                face_idx += n
                zone_idx += 1
        else:
            for b in range(nmf.blk_num):
                cur_face_total = blk_internal_face_num(nmf.blk[b].shape)
                face = np.empty((cur_face_total, 4), int)
                u, v, _ = nmf.blk[b].shape
                n = 0

                for j in range(1, v - 1):
                    for i in range(u - 1):
                        p = np.array([i, j])
                        crd_list = pnt_circle_h(p)
                        face[n][0] = nmf.calc_pnt_idx(b, crd_list[0])
                        face[n][1] = nmf.calc_pnt_idx(b, crd_list[1])
                        face[n][2] = nmf.calc_cell_idx(b, p, 1)
                        face[n][3] = nmf.calc_cell_idx(b, p, 4)
                        n += 1

                for i in range(1, u - 1):
                    for j in range(v - 1):
                        p = np.array([i, j])
                        crd_list = pnt_circle_v(p)
                        face[n][0] = nmf.calc_pnt_idx(b, crd_list[0])
                        face[n][1] = nmf.calc_pnt_idx(b, crd_list[1])
                        face[n][2] = nmf.calc_cell_idx(b, p, 2)
                        face[n][3] = nmf.calc_cell_idx(b, p, 1)
                        n += 1

                msh.add(XFComment("Blk-{} internal edges:".format(b)))
                msh.add(XFFace(zone_idx, face_idx, face_idx + n - 1, BCType.Interior, FaceType.Linear, face))
                face_idx += n
                zone_idx += 1

        '''Boundary Face'''
        face_type = FaceType.Quadrilateral if dim == 3 else FaceType.Linear
        for k, entry in enumerate(nmf.desc):
            adj_info = xf_calc_boundary_info(entry, nmf)
            bc_type = FluentMSH.NMF2MSH_BC_DICT[entry.type]
            msh.add(XFComment('Boundary-{}:'.format(k + 1)))
            msh.add(XFFace(zone_idx, face_idx, face_idx + n - 1, bc_type, face_type, adj_info))
            face_idx += n
            zone_idx += 1

        return msh
