import numpy as np
import math
from scipy.misc import comb
from scipy.linalg import solve
from src.nurbs.basis import equal, to_homogeneous, PntDist, find_span, Basis
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

        self.N = Basis(self.U, self.p)
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

    def w(self, u, d=0):
        ans = 0.0
        for i in range(0, self.n + 1):
            ans += self.N(i, self.p, u, d) * self.Pw[i][3]

        return ans

    def __call__(self, u, d=0):
        """
        计算曲线上的点
        :param u: 目标参数
        :param d: 求导次数
        :return: 曲线在u处的d阶导矢
        """

        Ad = np.zeros(3)
        nipd = np.zeros(self.n + 1)
        for i in range(0, self.n + 1):
            nipd[i] = self.N(i, self.p, u, d)
        for i in range(0, self.n + 1):
            for j in range(0, 3):
                Ad[j] += nipd[i] * self.Pw[i][j]

        ccw = np.zeros(3)
        i = 1
        while i <= d:
            ccw += comb(d, i, exact=True) * self.w(u, i) * self.__call__(u, d - i)
            i += 1

        wu = self.w(u)
        if equal(wu, 0.0):
            raise ValueError("Invalid weight settings!")

        return (Ad - ccw) / wu

    def to_iges(self, planar, periodic, norm, form=0):
        return IGES_Entity126(self.p, self.n, planar, (1 if self.isClosed else 0),
                              (1 if self.isPoly else 0), periodic,
                              self.U, self.weight, self.cpt,
                              self.U[0], self.U[-1], norm, form)

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
            Qw = calc_ctrl_pts(nU, p + 1, PPw, us)

            '''Update'''
            self.__init__(nU, Qw)
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
            alpha = 1.0 if i <= k - p else (0.0 if i >= k + 1 else (u - self.U[i]) / (self.U[i + p] - self.U[i]))
            if equal(alpha, 1.0):
                nPw[i] = self.Pw[i]
            elif equal(alpha, 0.0):
                nPw[i] = self.Pw[i - 1]
            else:
                nPw[i] = alpha * self.Pw[i] + (1 - alpha) * self.Pw[i - 1]

        '''Update'''
        self.__init__(nU, nPw)


class GlobalInterpolatedCrv(NURBS_Curve):
    def __init__(self, pts, p, method='centripetal'):
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
        dist[i] = PntDist(pts[i - 1], pts[i])

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

    n = len(param)
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
    N = Basis(U, p)
    cm = np.zeros((n + 1, n + 1))
    for k in range(0, n + 1):
        for i in range(0, n + 1):
            cm[k][i] = N(i, p, param[k])

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
