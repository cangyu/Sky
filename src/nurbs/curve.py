import numpy as np
import math
from numpy.linalg import norm
from scipy.linalg import solve
from scipy.interpolate import BSpline
from scipy.integrate import romberg
from src.nurbs.utility import *
from src.iges.iges_entity126 import IGES_Entity126
from src.iges.iges_entity110 import IGES_Entity110


class NURBS_Curve(object):
    def __init__(self, U, Pw):
        """
        NURBS曲线
        :param U:节点矢量 
        :param Pw: 带权控制点
        """

        self.n = len(Pw) - 1
        self.m = len(U) - 1
        self.p = self.m - self.n - 1

        self.U = np.copy(U)
        self.Pw = np.copy(Pw)
        self.spl = BSpline(self.U, self.Pw, self.p)

        self.weight = np.zeros(self.n + 1)
        self.cpt = np.zeros((self.n + 1, 3))
        self.isPoly = True
        self.isClosed = (self.__call__(self.U[0]) == self.__call__(self.U[-1])).all()

        for i in range(self.n + 1):
            self.weight[i] = self.Pw[i][3]
            if self.isPoly and (not equal(self.weight[i], 1.0)):
                self.isPoly = False
            for j in range(3):
                self.cpt[i][j] = Pw[i][j] / self.weight[i]

    def __call__(self, u, d=0):
        """
        计算曲线上的点
        :param u: 目标参数
        :param d: 求导次数
        :return: 曲线在u处的d阶导矢
        """

        pw = self.spl(u, d)
        return to_cartesian(pw)

    def length(self):
        """
        近似计算曲线长度
        """

        return romberg(lambda u: norm(self.__call__(u, 1)), 0.0, 1.0)

    def curvature(self, u):
        """
        计算在给定位置处曲率
        """

        p1 = self.__call__(u, 1)
        p2 = self.__call__(u, 2)

        dd = np.zeros(3)
        dd[0] = math.pow(p2[2] * p1[1] - p2[1] * p1[2], 2)
        dd[1] = math.pow(p2[0] * p1[2] - p2[2] * p1[0], 2)
        dd[2] = math.pow(p2[1] * p1[0] - p2[0] * p1[1], 2)
        dividend = math.sqrt(np.sum(dd))

        dv = np.zeros(3)
        dv[0] = math.pow(p1[0], 2)
        dv[1] = math.pow(p1[1], 2)
        dv[2] = math.pow(p1[2], 2)
        divisor = math.pow(np.sum(dv), 3 / 2)

        kappa = dividend / divisor
        return kappa

    def reset(self, U, Pw):
        self.n = len(Pw) - 1
        self.m = len(U) - 1
        self.p = self.m - self.n - 1
        self.U = np.copy(U)
        self.Pw = np.copy(Pw)
        self.spl = BSpline(self.U, self.Pw, self.p)
        self.weight = np.zeros(self.n + 1)
        self.cpt = np.zeros((self.n + 1, 3))
        self.isPoly = True
        self.isClosed = (self.__call__(self.U[0]) == self.__call__(self.U[-1])).all()
        for i in range(self.n + 1):
            self.weight[i] = self.Pw[i][3]
            if self.isPoly and (not equal(self.weight[i], 1.0)):
                self.isPoly = False
            for j in range(3):
                self.cpt[i][j] = Pw[i][j] / self.weight[i]

    def reverse(self):
        """
        曲线反向
        """

        nPw = self.Pw[::-1, :]
        nU = np.ones(self.m + 1) - self.U[::-1]
        self.reset(nU, nPw)

    def pan(self, delta):
        """
        平移
        :param delta: 偏移矢量
        :return: None
        """

        dv = np.zeros(3)
        array_smart_copy(delta, dv)
        nU = np.copy(self.U)
        nP = np.copy(self.cpt)
        nPw = np.empty(np.asarray(self.Pw))

        for i in range(self.n + 1):
            nP[i] += dv
            nPw[i] = to_homogeneous(nP[i], self.weight[i])

        self.reset(nU, nPw)

    def rotate(self):
        pass

    def to_iges(self, planar=0, periodic=0, norm_vector=np.zeros(3), form=0):
        return IGES_Entity126(self.p, self.n, planar, (1 if self.isClosed else 0), (1 if self.isPoly else 0), periodic,
                              self.U, self.weight, self.cpt, self.U[0], self.U[-1], norm_vector, form)

    def elevate(self, t: int):
        """
        将曲线升阶t次
        :param t: 升阶次数
        :return: None.
        """

        if t <= 0:
            return
        elif t == 1:
            '''New knot sequence'''
            val, cnt = np.unique(self.U, return_counts=True)
            s = len(val) - 2
            for i in range(0, len(cnt)):
                cnt[i] += 1

            nU = np.zeros(np.sum(cnt))
            k = 0
            for i in range(0, len(val)):
                for j in range(0, cnt[i]):
                    nU[k] = val[i]
                    k += 1

            '''Sample'''
            nn = self.n + 1 + s
            us = np.linspace(self.U[0], self.U[-1], nn + 1)

            PPw = np.zeros((nn + 1, 4))
            for i in range(0, len(us)):
                PPw[i] = to_homogeneous(self.__call__(us[i]))

            '''Solve'''
            Qw = calc_ctrl_pts(nU, self.p + 1, PPw, us)

            '''Update'''
            self.update(nU, Qw)
        else:
            self.elevate(t - 1)

    def insert_knot(self, u):
        """
        插入一个节点
        :param u: 待插入节点
        :return: None
        """

        '''Insert'''
        k = find_span(self.U, u)
        nU = np.insert(self.U, k, u)

        '''Calculate new CtrlPts'''
        n, dim = self.Pw.shape
        nPw = np.zeros((n + 1, dim))
        for i in range(0, n + 1):
            alpha = 1.0 if i <= k - self.p else (0.0 if i >= k + 1 else (u - self.U[i]) / (self.U[i + self.p] - self.U[i]))
            if equal(alpha, 1.0):
                nPw[i] = self.Pw[i]
            elif equal(alpha, 0.0):
                nPw[i] = self.Pw[i - 1]
            else:
                nPw[i] = alpha * self.Pw[i] + (1 - alpha) * self.Pw[i - 1]

        '''Update'''
        self.update(nU, nPw)


class GlobalInterpolatedCrv(NURBS_Curve):
    def __init__(self, pts, p=3, method='centripetal'):
        """
        构造一条p次非有理B样条曲线插值于pts
        :param pts: 待插值点序列
        :param p: 目标曲线次数
        :param method: 计算插值点参数的方法
        """

        n, dim = pts.shape
        n -= 1
        param = calc_pnt_param(pts, method)
        U = calc_knot_vector(param, p)
        P = calc_ctrl_pts(U, p, pts, param)

        Pw = np.zeros((n + 1, dim + 1))
        for i in range(0, n + 1):
            Pw[i] = to_homogeneous(P[i])

        super(GlobalInterpolatedCrv, self).__init__(U, Pw)


def calc_pnt_param(pts, method):
    """
    计算每个插值点所对应的参数。
    :param pts: 插值点坐标序列
    :param method: 参数计算方法
    :return: 插值点坐标序列对应的参数序列([0,1])
    """

    if method not in ['chord', 'centripetal']:
        raise ValueError("Invalid method parameter!")

    n = len(pts) - 1
    param = np.zeros(n + 1)
    param[n] = 1.0

    dist = np.zeros(n + 1)
    for i in range(1, n + 1):
        dist[i] = pnt_dist(pts[i - 1], pts[i])

    d = 0
    if method == 'chord':  # 弦长参数化
        for i in range(1, n + 1):
            d += dist[i]
    else:  # 向心参数化，数据点急转弯变化时效果好
        for i in range(1, n + 1):
            dist[i] = math.sqrt(dist[i])
            d += dist[i]

    for i in range(1, n):
        param[i] = param[i - 1] + dist[i] / d

    return param


def calc_knot_vector(param, p):
    """
    取平均值方法计算节点
    :param param: 插值点序列对应的参数序列
    :param p: 目标曲线次数
    :return: 目标曲线节点矢量([0,1])
    """

    n = len(param) - 1
    m = n + p + 1
    knots = np.zeros(m + 1)

    '''Tail'''
    for i in range(0, p + 1):
        knots[m - i] = 1.0

    '''Prepare'''
    acc = 0.0
    for i in range(0, p):
        acc += param[i]

    '''Iterate'''
    for j in range(1, n - p + 1):
        acc -= param[j - 1]
        acc += param[p - 1 + j]
        knots[p + j] = acc / p

    return knots


def calc_ctrl_pts(U, p, pts, param):
    """
    求解线性方程组得到控制点
    :param U: 节点矢量
    :param p: 目标曲线次数
    :param pts: 插值点序列
    :param param: 插值点所对应参数
    :return: 控制点序列
    """
    n, dim = pts.shape
    n -= 1

    ctrl_pts = np.zeros((n + 1, dim))

    '''Coefficient Matrix'''
    cm = np.zeros((n + 1, n + 1))
    for k in range(0, n + 1):
        cm[k] = all_basis_val(param[k], p, U)

    '''Solve'''
    Q = np.zeros((dim, n + 1))
    P = np.zeros((dim, n + 1))

    for i in range(0, dim):
        for j in range(0, n + 1):
            Q[i][j] = pts[j][i]

    for i in range(0, dim):
        P[i] = solve(cm, Q[i])

    for i in range(0, n + 1):
        for j in range(0, dim):
            ctrl_pts[i][j] = P[j][i]

    return ctrl_pts


class Line(NURBS_Curve):
    def __init__(self, a, b):
        """
        两点间直线段
        :param a: 起始点坐标
        :param b: 终点坐标
        """

        U = np.array([0, 0, 1, 1])
        Pw = np.array([[0, 0, 0, 1], [0, 0, 0, 1]], float)
        array_smart_copy(a, Pw[0])
        array_smart_copy(b, Pw[1])

        super(Line, self).__init__(U, Pw)

    def length(self):
        return pnt_dist(to_cartesian(self.Pw[0]), to_cartesian(self.Pw[-1]))

    def curvature(self, u):
        return 0.0

    def to_iges(self):
        return IGES_Entity110(to_cartesian(self.Pw[0]), to_cartesian(self.Pw[-1]))


class Arc(NURBS_Curve):
    def __init__(self, r, theta):
        """
        XY平面内简化圆弧，以原点为圆心，起始点为(r,0),法向量为(0,0,1)
        :param r: 半径
        :param theta: 圆心角, will be fitted into range (0,360]
        """

        while theta <= 0:
            theta += 360
        while theta > 360:
            theta -= 360

        self.radius = r
        self.theta = theta

        narcs = int(np.ceil(theta / 90))
        theta = np.deg2rad(theta)
        dtheta = theta / narcs
        w1 = np.cos(dtheta / 2)
        dknot = 1.0 / narcs

        m = 2 * narcs + 3  # 最后一个节点下标
        n = 2 * narcs  # 最后一个控制点下标

        U = np.zeros(m + 1)
        P = np.zeros((n + 1, 3))
        Pw = np.zeros((n + 1, 4))

        '''Knot Vector'''
        U[-1] = U[-2] = U[-3] = 1.0
        for i in range(1, narcs):
            cur_index = 1 + 2 * i
            U[cur_index] = U[cur_index + 1] = i * dknot

        '''Control Points'''
        self.sp = P[0] = np.array([r, 0, 0], float)
        Pw[0] = to_homogeneous(P[0], 1.0)
        T0 = np.array([0.0, 1.0, 0.0])

        index = 0
        angle = 0.0
        for i in range(1, narcs + 1):
            angle += dtheta
            P[index + 2] = np.array([r * np.cos(angle), r * np.sin(angle), 0.0])
            Pw[index + 2] = to_homogeneous(P[index + 2], 1.0)
            T2 = np.array([-np.sin(angle), np.cos(angle), 0.0])
            P[index + 1] = line_intersection(P[index], T0, P[index + 2], T2)
            Pw[index + 1] = to_homogeneous(P[index + 1], w1)
            index += 2
            if i < narcs:
                T0 = T2

        self.ep = P[-1]
        super(Arc, self).__init__(U, Pw)

    def curvature(self, u):
        return 1.0 / self.radius

    def length(self):
        return self.radius * np.deg2rad(self.theta)

    @classmethod
    def from_2pnt(cls, start_pnt, end_pnt, theta, norm_vector):
        """
        空间圆弧
        :param start_pnt: 起始点坐标
        :param end_pnt: 终止点坐标
        :param theta: 圆心角
        :param norm_vector: 所在平面的法向量，按右手定则， 四指依次扫过start_pnt和end_pnt
        """

        '''Basic Variables'''
        sp = np.copy(start_pnt)
        ep = np.copy(end_pnt)
        theta = np.deg2rad(theta)
        radius = 0.5 * pnt_dist(sp, ep) / np.sin(theta / 2)
        w = radius * np.cos(theta / 2)
        cdir = np.cross(norm_vector, ep - sp)
        cdir /= norm(cdir)
        origin = 0.5 * (sp + ep) + cdir * w

        arc = cls(radius, np.rad2deg(theta))
        arc.pan(origin)
        arc.rotate()
        return arc
