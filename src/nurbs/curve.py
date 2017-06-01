import numpy as np
import math
from scipy.linalg import solve
from scipy.interpolate import BSpline
from src.nurbs.utility import to_cartesian, to_homogeneous, pnt_dist, equal, all_basis_val
from src.iges.iges_entity126 import IGES_Entity126


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

        for i in range(0, self.n + 1):
            self.weight[i] = self.Pw[i][3]
            if self.isPoly and (not equal(self.weight[i], 1.0)):
                self.isPoly = False
            for j in range(0, 3):
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

    def to_iges(self, planar, periodic, norm, form=0):
        return IGES_Entity126(self.p, self.n, planar, (1 if self.isClosed else 0),
                              (1 if self.isPoly else 0), periodic,
                              self.U, self.weight, self.cpt,
                              self.U[0], self.U[-1], norm, form)


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
