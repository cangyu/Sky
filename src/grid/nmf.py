#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

"""
Implementation of NASA's Neutral Map File representation.

Ref:
https://geolab.larc.nasa.gov/Volume/Doc/nmf.htm
"""


def blk_node_num(shape):
    """
    Calculate the num of nodes in a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of nodes in the block.
    :rtype: int
    """

    if len(shape) == 3:
        u, v, _ = shape
        return u * v
    elif len(shape) == 4:
        u, v, w, _ = shape
        return u * v * w
    else:
        raise ValueError('Invalid shape.')


def blk_face_num(shape):
    """
    Calculate the num of faces in a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of faces in the block.
    :rtype: int
    """

    if len(shape) == 3:
        u, v, _ = shape
        return u * (v - 1) + v * (u - 1)
    elif len(shape) == 4:
        u, v, w, _ = shape
        return u * (v - 1) * (w - 1) + v * (w - 1) * (u - 1) + w * (u - 1) * (v - 1)
    else:
        raise ValueError('Invalid shape.')


def blk_cell_num(shape):
    """
    Calculate the num of cells in a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of cells in the block.
    :rtype: int
    """

    if len(shape) == 3:
        u, v, _ = shape
        return (u - 1) * (v - 1)
    elif len(shape) == 4:
        u, v, w, _ = shape
        return (u - 1) * (v - 1) * (w - 1)
    else:
        raise ValueError('Invalid shape.')


def blk_internal_node_num(shape):
    """
    Calculate the num of nodes inside a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of cells inside the block.
    :rtype: int
    """

    if len(shape) == 3:
        u, v, _ = shape
        return (u - 2) * (v - 2)
    elif len(shape) == 4:
        u, v, w, _ = shape
        return (u - 2) * (v - 2) * (w - 2)
    else:
        raise ValueError('Invalid shape.')


def blk_internal_face_num(shape):
    """
    Calculate the num of faces inside a block.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: Num of faces inside the block.
    :rtype: int
    """

    if len(shape) == 3:
        u, v, _ = shape
        return (u - 2) * (v - 1) + (u - 1) * (v - 2)
    elif len(shape) == 4:
        u, v, w, _ = shape
        return (u - 1) * (v - 1) * (w - 2) + (u - 1) * (v - 2) * (w - 1) + (u - 2) * (v - 1) * (w - 1)
    else:
        raise ValueError('Invalid shape.')


def on_boundary(pnt, shape):
    """
    Judge if the point is on the boundary of the block.
    :param pnt: Coordinate of the point.
    :param shape: Shape of the block(Dimension of point is also included for convenience).
    :return: If the point is on the boundary.
    :rtype: bool
    """

    assert len(pnt) == len(shape) - 1

    if len(shape) == 3:
        i, j = pnt
        u, v, _ = shape
        return i == 0 or i == u - 1 or j == 0 or j == v - 1
    elif len(shape) == 4:
        i, j, k = pnt
        u, v, w, _ = shape
        return i == 0 or i == u - 1 or j == 0 or j == v - 1 or k == 0 or j == w - 1
    else:
        raise ValueError('Invalid shape.')


def blk_node_idx(pnt, shape):
    """
    The node index in a single block(2D or 3D).
    :param pnt: Coordinate of the point.
    :param shape: Shape of the block.
    :return: Index with (I > J > K) preference, starting from 0.
    :rtype: int
    """

    assert len(pnt) == len(shape) - 1

    if len(shape) == 3:
        u, v, _ = shape
        i, j = pnt
        return j * u + i
    elif len(shape) == 4:
        u, v, w, _ = shape
        i, j, k = pnt
        return k * u * v + j * u + i
    else:
        raise ValueError('Invalid shape.')


def blk_internal_node_idx(pnt, shape):
    """
    The node index within a block, counting from the internal.
    :param pnt: Coordinate of the point.
    :param shape: Shape of the block.
    :return: Internal index of the point with (I > J > K) preference, starting from 0.
    :rtype: int
    """

    assert len(pnt) == len(shape) - 1

    if len(shape) == 3:
        u, v, _ = shape
        i, j = pnt
        return blk_node_idx((i - 1, j - 1), (u - 2, v - 2))
    elif len(shape) == 4:
        u, v, w, _ = shape
        i, j, k = pnt
        return blk_node_idx((i - 1, j - 1, k - 1), (u - 2, v - 2, w - 2))
    else:
        raise ValueError('Invalid shape.')


def boundary_pnt_caste(shape):
    """
    Calculate the caste of the boundary points on a single block.
    Used to serialize all the boundary points.
    :param shape: The shape of the block.
    :return: Caste of the boundary points.
    """

    u, v, w, _ = shape
    ret = np.empty(6, int)
    ret[0] = ret[1] = w * v
    ret[2] = ret[3] = (u - 2) * w
    ret[4] = ret[5] = (u - 2) * (v - 2)
    for i in range(1, 6):
        ret[i] += ret[i - 1]
    return ret


def blk_cell_idx(pnt, shape):
    """
    Calculate the cell index in a block through its corner point.
    :param pnt: Coordinate of the corner point.
    :param shape: The shape of the block.
    :return: Take the corner point as origin, inherit the positive directions from block,
             calculate the cell index at the 1st quadrant(Starting from 0).
    :rtype: int
    """

    assert len(pnt) == len(shape) - 1

    if len(shape) == 4:
        u, v, w, _ = shape
        i, j, k = pnt
        assert 0 <= i < u - 1 and 0 <= j < v - 1 and 0 <= k < w - 1
        return (u - 1) * (v - 1) * k + (u - 1) * j + i
    elif len(shape) == 3:
        u, v, _ = shape
        i, j = pnt
        assert 0 <= i < u - 1 and 0 <= j < v - 1
        return j * (u - 1) + i
    else:
        raise ValueError('Invalid shape.')


def blk_cell_idx_quadrant(pnt, shape, q):
    """
    Calculate the cell index in a block through its corner point and quadrant number.
    :param pnt: Coordinate of the corner point.
    :param shape: The shape of the block.
    :param q: Quadrant number.
    :type q: int
    :return: Take the corner point as origin, inherit the positive directions from block,
             calculate the index of cell at specified quadrant(Starting from 0).
    :rtype: int
    """

    assert len(pnt) == len(shape) - 1

    if len(shape) == 4:
        i, j, k = pnt
        if q == 1:
            return blk_cell_idx(pnt, shape)
        elif q == 2:
            return blk_cell_idx((i - 1, j, k), shape)
        elif q == 3:
            return blk_cell_idx((i - 1, j - 1, k), shape)
        elif q == 4:
            return blk_cell_idx((i, j - 1, k), shape)
        elif q == 5:
            return blk_cell_idx((i, j, k - 1), shape)
        elif q == 6:
            return blk_cell_idx((i - 1, j, k - 1), shape)
        elif q == 7:
            return blk_cell_idx((i - 1, j - 1, k - 1), shape)
        elif q == 8:
            return blk_cell_idx((i, j - 1, k - 1), shape)
        else:
            raise ValueError('Invalid quadrant number.')
    elif len(shape) == 3:
        assert 1 <= q <= 4
        i, j = pnt
        if q == 1:
            return blk_cell_idx(pnt, shape)
        elif q == 2:
            return blk_cell_idx((i - 1, j), shape)
        elif q == 3:
            return blk_cell_idx((i - 1, j - 1), shape)
        elif q == 4:
            return blk_cell_idx((i, j - 1), shape)
        else:
            raise ValueError('Invalid quadrant number.')
    else:
        raise ValueError('Invalid shape.')


def face_invariant_idx(f, shape):
    """
    Get the invariant index on certain face according to NMF convention.
    :param f: Face number.
    :type f: int
    :param shape: Shape of the block.
    :return: Invariant index for face 'f'.
    """

    u, v, w, _ = shape
    if f == 1 or f == 3 or f == 5:
        return 1
    elif f == 2:
        return w
    elif f == 4:
        return u
    elif f == 6:
        return v
    else:
        raise ValueError('Invalid face number.')


class NMFEntry(object):
    ONE_TO_ONE = 'ONE_TO_ONE'
    NMF_LITERAL = '# {:<18}{:^8}{:^2}{:^8}{:^8}{:^8}{:^8}{:^8}{:^2}{:^8}{:^8}{:^8}{:^8}{:^6}'.format('Type', 'B1', 'F1', 'S1', 'E1', 'S2', 'E2', 'B2', 'F2', 'S1', 'E1', 'S2', 'E2', 'Swap')
    IDX_MAP = {1: (2, 0, 1), 2: (2, 0, 1), 3: (0, 1, 2), 4: (0, 1, 2), 5: (1, 2, 0), 6: (1, 2, 0)}
    # CornerIdx in b2, Pri_Rev, Sec_Rev, Swp
    FACE_MAPPING = [[0, 1, 2, 3, False, False, False],
                    [0, 3, 2, 1, False, False, True],
                    [1, 0, 3, 2, True, False, False],
                    [3, 0, 1, 2, True, False, True],
                    [2, 3, 0, 1, True, True, False],
                    [2, 1, 0, 3, True, True, True],
                    [3, 2, 1, 0, False, True, False],
                    [1, 2, 3, 0, False, True, True]]

    def __init__(self, *args, **kwargs):
        """
        Entry used in Neutral Map File to indicate the topological features of the mesh.
        :param args: Parameters.
        :param kwargs: Options.
        """

        assert len(args) in (7, 13)

        self.Type = args[0]  # The type of feature (topological or boundary condition) to be defined and positioned within the mesh.
        self.B1 = args[1]  # The number of the first block(Starting from 1).
        self.F1 = args[2]  # The face number for the first block(From 1 to 6).
        self.B1PriStart = args[3]  # The starting index in the primary coordinate direction for the face in the 1st block.
        self.B1PriEnd = args[4]  # The ending index in the primary coordinate direction for the face in the 1st block.
        self.B1SecStart = args[5]  # The starting index in the secondary coordinate direction for the face in the 1st block.
        self.B1SecEnd = args[6]  # The ending index in the secondary coordinate direction for the face in the 1st block.
        assert self.B1 > 0 and 1 <= self.F1 <= 6
        self.B1Shape = np.zeros(4, int)

        if len(args) == 13:
            self.B2 = args[7]  # The number of the second block(Starting from 1).
            self.F2 = args[8]  # The face number for the second block(From 1 to 6).
            self.B2PriStart = args[9]  # The starting index in the primary coordinate direction for the face in the 2nd block.
            self.B2PriEnd = args[10]  # The ending index in the primary coordinate direction for the face in the 2nd block.
            self.B2SecStart = args[11]  # The starting index in the secondary coordinate direction for the face in the 2nd block.
            self.B2SecEnd = args[12]  # The ending index in the secondary coordinate direction for the face in the 2nd block.
            assert self.B2 > 0 and 1 <= self.F2 <= 6
            self.B2Shape = np.zeros(4, int)

            '''
            Orientation flag(Specified only for Type==ONE_TO_ONE).
                False - The primary directions of the two faces are aligned
                True - Otherwise
            At this stage, this setting may be not correct, will be automatically configured later.
            '''
            self.Swap = kwargs['swap'] if 'swap' in kwargs else False

            '''
            Even directions are aligned, they may be opposite.
            At this stage, these settings may be not correct, will be automatically configured later.
            '''
            self.PriReverse = False
            self.SecReverse = False

        '''Indicate relative position'''
        self.OnLeft = True

        '''Points back to parent'''
        self.NMF = None

    '''
    Followings are some getters and setters
    for easy use.
    '''

    @property
    def shape_of_blk1(self):
        return self.B1Shape

    @shape_of_blk1.setter
    def shape_of_blk1(self, shape):
        assert len(shape) == len(self.B1Shape)
        self.B1Shape = np.copy(shape)

    @property
    def shape_of_blk2(self):
        return self.B2Shape

    @shape_of_blk2.setter
    def shape_of_blk2(self, shape):
        assert len(shape) == len(self.B2Shape)
        self.B2Shape = np.copy(shape)

    def shape_of_blk(self, b=1):
        if b == 1:
            return self.shape_of_blk1
        elif b == 2:
            return self.shape_of_blk2
        else:
            raise ValueError('invalid block indication')

    '''
    Followings class function are used
    as constructors.
    '''

    @classmethod
    def single(cls, tp, b1, f1, s1, e1, s2, e2):
        return cls(tp, b1, f1, s1, e1, s2, e2)

    @classmethod
    def one2one(cls, *args):
        if len(args) == 13:
            b1, f1, b1s1, b1e1, b1s2, b1e2, b2, f2, b2s1, b2e1, b2s2, b2e2, swp = args
            return cls(NMFEntry.ONE_TO_ONE, b1, f1, b1s1, b1e1, b1s2, b1e2, b2, f2, b2s1, b2e1, b2s2, b2e2, swap=swp)
        elif len(args) == 12:
            b1, f1, b1s1, b1e1, b1s2, b1e2, b2, f2, b2s1, b2e1, b2s2, b2e2 = args
            return cls(NMFEntry.ONE_TO_ONE, b1, f1, b1s1, b1e1, b1s2, b1e2, b2, f2, b2s1, b2e1, b2s2, b2e2)
        else:
            raise ValueError('invalid arguments')

    '''
    Following utilities are used for
    ASCII output and representation.
    '''

    def __repr__(self):
        ret = '{:<20}'.format(self.Type)
        ret += '{:^8}'.format(self.B1)
        ret += '{:^2}'.format(self.F1)
        ret += '{:^8}'.format(self.B1PriStart)
        ret += '{:^8}'.format(self.B1PriEnd)
        ret += '{:^8}'.format(self.B1SecStart)
        ret += '{:^8}'.format(self.B1SecEnd)
        if self.Type == NMFEntry.ONE_TO_ONE:
            ret += '{:^8}'.format(self.B2)
            ret += '{:^2}'.format(self.F2)
            ret += '{:^8}'.format(self.B2PriStart)
            ret += '{:^8}'.format(self.B2PriEnd)
            ret += '{:^8}'.format(self.B2SecStart)
            ret += '{:^8}'.format(self.B2SecEnd)
            ret += '{:>6}'.format('TRUE' if self.Swap else 'FALSE')
        return ret

    def write(self, f_out):
        f_out.write(self.__repr__())

    '''
    Following routines and properties
    are used for basic counting.
    '''

    def pri_node_num(self, b=1):
        if b == 1:
            return self.B1PriEnd - self.B1PriStart + 1
        elif b == 2:
            return self.B2PriEnd - self.B2PriStart + 1
        else:
            raise ValueError('invalid block indication')

    def sec_node_num(self, b=1):
        if b == 1:
            return self.B1SecEnd - self.B1SecStart + 1
        elif b == 2:
            return self.B2SecEnd - self.B2SecStart + 1
        else:
            raise ValueError('invalid block indication')

    @property
    def node_num(self):
        t1 = self.pri_node_num()
        t2 = self.sec_node_num()
        return t1 if t2 == 0 else t1 * t2

    @property
    def face_num(self):
        t1 = self.pri_node_num() - 1
        t2 = self.sec_node_num() - 1
        return t1 if t2 == 0 else t1 * t2

    '''
    Followings are utilities to calculate
    the invariant index of each side.
    Note that the result starting from 1.
    '''

    @property
    def f1_inv_idx(self):
        return face_invariant_idx(self.F1, self.B1Shape)

    @property
    def f2_inv_idx(self):
        return face_invariant_idx(self.F2, self.B2Shape)

    def invariant_idx(self, b=1):
        if b == 1:
            return self.f1_inv_idx
        elif b == 2:
            return self.f2_inv_idx
        else:
            raise ValueError('Invalid block indication.')

    '''
    Followings are routines to calculate
    the physical index corresponding to
    each logical axis.
    '''

    def b1_pri_idx(self, x1):
        return self.B1PriStart + x1

    def b1_sec_idx(self, x2):
        return self.B1SecStart + x2

    def b2_pri_idx(self, x1):
        return self.B2PriStart + x1

    def b2_sec_idx(self, x2):
        return self.B2SecStart + x2

    def pri_idx(self, x1, b=1):
        if b == 1:
            return self.b1_pri_idx(x1)
        elif b == 2:
            return self.b2_pri_idx(x1)
        else:
            raise ValueError('invalid block indication')

    def sec_idx(self, x2, b=1):
        if b == 1:
            return self.b1_sec_idx(x2)
        elif b == 2:
            return self.b2_sec_idx(x2)
        else:
            raise ValueError('invalid block indication')

    '''
    Followings are utilities for converting
    logical index into real index and vice-versa.
    Note that the primary axis on each face may be
    swapped, so the input is block dependent.
    '''

    def logic2real(self, lp, b=1):
        ret = np.empty(3, int)
        x1, x2 = lp
        if b == 1:
            mp = NMFEntry.IDX_MAP[self.F1]
            ret[mp[0]] = self.f1_inv_idx - 1
            ret[mp[1]] = self.b1_pri_idx(x1) - 1
            ret[mp[2]] = self.b1_sec_idx(x2) - 1
        elif b == 2:
            mp = NMFEntry.IDX_MAP[self.F2]
            ret[mp[0]] = self.f2_inv_idx - 1
            ret[mp[1]] = self.b2_pri_idx(x1) - 1
            ret[mp[2]] = self.b2_sec_idx(x2) - 1
        else:
            raise ValueError('Invalid block indication.')

        return ret

    def real2logic(self, rp, b=1):
        assert b in (1, 2)

        f = self.F1 if b == 1 else self.F2
        mp = NMFEntry.IDX_MAP[f]
        x = [1 + rp[mp[i]] for i in range(3)]

        if b == 1:
            x[1] -= self.B1PriStart
            x[2] -= self.B1SecStart
        else:
            x[1] -= self.B2PriStart
            x[2] -= self.B2SecStart

        assert x[0] == self.invariant_idx(b)
        return np.array([x[1], x[2]])

    '''
    Following routines are used to determine
    the direction mapping relations between
    the 2 faces automatically.
    '''

    def corners(self, b=1):
        ret = np.empty((4, 3), int)

        if b == 1:
            pnt = [(self.B1PriStart, self.B1SecStart),
                   (self.B1PriEnd, self.B1SecStart),
                   (self.B1PriEnd, self.B1SecEnd),
                   (self.B1PriStart, self.B1SecEnd)]
        elif b == 2:
            pnt = [(self.B2PriStart, self.B2SecStart),
                   (self.B2PriEnd, self.B2SecStart),
                   (self.B2PriEnd, self.B2SecEnd),
                   (self.B2PriStart, self.B2SecEnd)]
        else:
            raise ValueError('invalid block indication')

        mp = NMFEntry.IDX_MAP[self.F1 if b == 1 else self.F2]
        for i in range(4):
            ret[i][mp[0]] = self.invariant_idx(b) - 1
            ret[i][mp[1]] = pnt[i][0] - 1
            ret[i][mp[2]] = pnt[i][1] - 1

        return ret

    def direction_auto_mapping(self):
        c1 = self.corners(1)
        c2 = self.corners(2)

        ok = False
        for dt in NMFEntry.FACE_MAPPING:
            all_match = True
            for k in range(4):
                i1, j1, k1 = c1[k]
                i2, j2, k2 = c2[dt[k]]
                p1 = self.NMF.blk[self.B1 - 1][i1][j1][k1]
                p2 = self.NMF.blk[self.B2 - 1][i2][j2][k2]
                if not np.array_equal(p1, p2):
                    all_match = False
                    break
            if all_match:
                self.PriReverse = dt[4]
                self.SecReverse = dt[5]
                self.Swap = dt[6]
                ok = True
                break

        if not ok:
            raise ValueError('2 faces not match')

        '''Checking'''
        if self.Swap:
            assert self.pri_node_num(1) == self.sec_node_num(2) and self.sec_node_num(1) == self.pri_node_num(2)
        else:
            assert self.pri_node_num(1) == self.pri_node_num(2) and self.sec_node_num(1) == self.sec_node_num(2)

    '''
    Following routine is used for finding the
    counterpart logic index on the opposite side.
    '''

    def counterpart(self, lp, b):
        """
        Calculate corresponding logic coordinate in the opposite block.
        :param lp: Logic point.
        :param b: For block indication: 1-Blk1, 2-Blk2.
        :type b: int
        :return: Counterpart logic coordinate.
        """

        x1, y1 = lp
        x2 = self.pri_node_num(b) - 1 - x1 if self.PriReverse else x1
        y2 = self.sec_node_num(b) - 1 - y1 if self.SecReverse else y1
        if self.Swap:
            x2, y2 = y2, x2

        return np.array([x2, y2])


class NeutralMapFile(object):
    FACE_INVARIANT = {1: 2, 2: 2, 3: 0, 4: 0, 5: 1, 6: 1}
    CELL_QUADRANT_ON_FACE = {1: 1, 2: 5, 3: 1, 4: 2, 5: 1, 6: 4}

    def __init__(self, str_grid):
        self.blk = str_grid
        self.desc = []

        self.internal_pnt_num = np.array([blk_internal_node_num(self.blk[i].shape) for i in range(self.blk_num)])
        self.boundary_pnt_num = np.array([blk_node_num(self.blk[i].shape) for i in range(self.blk_num)])
        self.boundary_pnt_num -= self.internal_pnt_num
        for i in range(1, self.blk_num):
            self.internal_pnt_num[i] += self.internal_pnt_num[i - 1]
            self.boundary_pnt_num[i] += self.boundary_pnt_num[i - 1]

        self.boundary_pnt_idx = [np.zeros(self.boundary_pnt_num[i], int) for i in range(self.blk_num)]
        self.boundary_pnt_cnt = 0
        self.boundary_pnt_map = {}

        self.boundary_pnt_caste = np.array([boundary_pnt_caste(self.blk[i].shape) for i in range(self.blk_num)])

        self.cell_start = np.array([blk_cell_num(self.blk[i].shape) for i in range(self.blk_num)])
        for i in range(1, self.blk_num):
            self.cell_start[i] += self.cell_start[i - 1]

        self.internal_face_num = np.array([blk_internal_face_num(self.blk[i].shape) for i in range(self.blk_num)])
        for i in range(1, self.blk_num):
            self.internal_face_num[i] += self.internal_face_num[i - 1]

    def add(self, entry):
        """
        Add topological or boundary info of the grid.
        :param entry: Description of the topological info.
        :type entry: NMFEntry
        :return: None.
        """

        self.desc.append(entry)

    def save(self, fn):
        """
        Save the boundary mapping info into file.
        :param fn: File name.
        :type fn: str
        :return: None.
        """

        f_out = open(fn, 'w')
        f_out.write('# =========================== Neutral Map File generated by the Sky software =================================')
        f_out.write('# ============================================================================================================')
        f_out.write('# Block#    IDIM    JDIM    KDIM')
        f_out.write('# ------------------------------------------------------------------------------------------------------------')
        f_out.write('{:>8}'.format(self.blk_num))
        for i in range(self.blk_num):
            u, v, w, _ = self.blk[i].shape
            f_out.write('{:>8}{:>8}{:>8}{:>8}'.format(i + 1, u, v, w))
        f_out.write('# ============================================================================================================')
        f_out.write(NMFEntry.NMF_LITERAL)
        f_out.write('# ------------------------------------------------------------------------------------------------------------')
        for i in range(len(self.desc)):
            self.desc[i].write(f_out)
        f_out.close()

    @property
    def grid(self):
        return self.blk

    @property
    def dim(self):
        return len(self.blk[0].shape) - 1

    @property
    def blk_num(self):
        return len(self.blk)

    @property
    def cell_num(self):
        return self.cell_start[-1]

    @property
    def pnt_num(self):
        return self.internal_pnt_num[-1] + self.boundary_pnt_cnt

    @property
    def face_num(self):
        return self.internal_face_num[-1] + sum([e.face_num for e in self.desc])

    def compute_topology(self):
        """
        Indexing all nodes without duplication.
        :return: None.
        """

        n = self.boundary_pnt_num[-1]

        '''Nodes within a boundary'''
        for entry in self.desc:
            for x1 in range(1, entry.pri_node_num - 1):
                for x2 in range(1, entry.sec_node_num - 1):
                    p1 = self.calc_real_pnt(entry, (x1, x2), 1)
                    t1 = self.calc_boundary_pnt_seq(entry.B1, p1)
                    self.boundary_pnt_idx[entry.B1][t1] = n
                    if entry.B2 != 0:
                        p2 = self.calc_real_pnt(entry, (x1, x2), 2)
                        t2 = self.calc_boundary_pnt_seq(entry.B2, p2)
                        self.boundary_pnt_idx[entry.B2][t2] = n
                    n += 1

        '''Nodes on the boundary'''
        for entry in self.desc:
            lp = entry.all_boundary_logic_pnt
            for logic_pnt in lp:
                # For each point on the boundary, if it has been indexed, just skip;
                # if not, find all pnt related to this and index them all.
                real_pnt = self.calc_real_pnt(entry, logic_pnt, 1)
                seq = self.calc_boundary_pnt_seq(entry.B1, real_pnt)
                if self.boundary_pnt_idx[entry.B1][seq] == 0:
                    t = self.find_all_occurrence(entry.B1, real_pnt)
                    for item in t:
                        b, crd = item
                        seq = self.calc_boundary_pnt_seq(b, crd)
                        self.boundary_pnt_idx[b][seq] = n
                    n += 1

        self.boundary_pnt_cnt = n - self.boundary_pnt_num[-1]

    def find_all_occurrence(self, b, p):
        ret = [(b, p)]
        t = 0
        while t < len(ret):
            cb, cp = ret[t]
            face = boundary_pnt_face(cp, self.blk[cb].shape)
            for f in face:
                if adj_desc[cb][f - 1][1] != 0:
                    b2, f2, swap = adj_desc[cb][f - 1]
                    cp_crd = get_counterpart_pnt_coord(f, b2, f2, crd, swap)
                    p2 = shell_pnt_idx_from_coord(b2, cp_crd)
                    ca = (b2, p2)
                    if ca not in ret:
                        ret.append(ca)
            t += 1

        return ret

    @property
    def all_pnt(self):
        ret = np.empty((self.pnt_num, 3), float)
        n = 0
        for b in range(self.blk_num):
            u, v, w, _ = self.blk[b].shape
            for k in range(1, w - 1):
                for j in range(1, v - 1):
                    for i in range(1, u - 1):
                        ret[n] = self.blk[b][i][j][k]
                        n += 1

        for t in range(self.boundary_pnt_cnt):
            b, crd = self.boundary_pnt_map[n][0]
            i, j, k = crd
            ret[n] = self.blk[b][i][j][k]
            n += 1

        return ret

    def calc_boundary_pnt_seq(self, b, pnt):
        u, v, w, _ = self.blk[b].shape
        i, j, k = pnt

        if i == 0:
            return j + k * v
        elif i == u - 1:
            return self.boundary_pnt_caste[b][0] + j + k * v
        else:
            if j == 0:
                return self.boundary_pnt_caste[b][1] + k + (i - 1) * w
            elif j == v - 1:
                return self.boundary_pnt_caste[b][2] + k + (i - 1) * w
            else:
                if k == 0:
                    return self.boundary_pnt_caste[b][3] + (i - 1) + (j - 1) * (u - 2)
                elif k == w - 1:
                    return self.boundary_pnt_caste[b][4] + (i - 1) + (j - 1) * (u - 2)
                else:
                    raise ValueError("Not a boundary pnt.")

    def calc_pnt_idx(self, b, pnt):
        """
        Given a point and its blk idx, calculate the global index of that point.
        :param b: Block index.
        :type b: int
        :param pnt: Coordinate of the point.
        :return: The global index of the point.
        :rtype: int
        """

        shape = self.blk[b].shape
        if on_boundary(pnt, shape):
            t = self.calc_boundary_pnt_seq(b, pnt)
            return self.boundary_pnt_idx[b][t]
        else:
            base = 0 if b == 0 else self.internal_pnt_num[b - 1]
            off = blk_internal_node_idx(pnt, shape)
            return base + off

    def calc_boundary_pnt(self, b, seq):
        u, v, w, _ = self.blk[b].shape
        if seq < self.boundary_pnt_caste[b][0]:
            i = 0
            j = seq % v
            k = seq // v
        elif seq < self.boundary_pnt_caste[b][1]:
            cp = seq - self.boundary_pnt_caste[b][0]
            i = u - 1
            j = cp % v
            k = cp // v
        elif seq < self.boundary_pnt_caste[b][2]:
            cp = seq - self.boundary_pnt_caste[b][1]
            i = 1 + cp // w
            j = 0
            k = cp % w
        elif seq < self.boundary_pnt_caste[b][3]:
            cp = seq - self.boundary_pnt_caste[b][2]
            i = 1 + cp // w
            j = v - 1
            k = cp % w
        elif seq < self.boundary_pnt_caste[b][4]:
            cp = seq - self.boundary_pnt_caste[b][3]
            i = cp % (u - 2) + 1
            j = cp // (u - 2) + 1
            k = 0
        elif seq < self.boundary_pnt_caste[b][5]:
            cp = seq - self.boundary_pnt_caste[b][4]
            i = cp % (u - 2) + 1
            j = cp // (u - 2) + 1
            k = w - 1
        else:
            raise ValueError("Invalid pnt index.")

        return np.array([i, j, k])

    def calc_cell_idx(self, b, pnt, quadrant):
        base = 0 if b == 0 else self.cell_start[b - 1]
        off = blk_cell_idx_quadrant(pnt, self.blk[b].shape, quadrant)
        return base + off

    def cell_between_face(self, b, rp1, rp2):
        pass
