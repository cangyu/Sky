from numpy.linalg import norm
from scipy.linalg import solve
from scipy.interpolate import BSpline
from scipy.integrate import romberg
from scipy.misc import comb
from src.nurbs.utility import *
from src.transform.dcm import DCM, normalize
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

        sp = self.__call__(self.U[0])
        ep = self.__call__(self.U[-1])
        self.isClosed = (sp == ep).all()

        for i in range(self.n + 1):
            self.weight[i] = self.Pw[i][3]
            if self.isPoly and (not equal(self.weight[i], 1.0)):
                self.isPoly = False
            for j in range(3):
                self.cpt[i][j] = Pw[i][j] / self.weight[i]

    def __call__(self, u, d=0, return_cartesian=True):
        """
        计算曲线上的点
        :param u: 目标参数
        :param d: 求导次数
        :return: 曲线在u处的d阶导矢
        """

        pw = self.spl(u, d)
        return to_cartesian(pw) if return_cartesian else pw

    def length(self):
        """
        近似计算曲线长度
        """

        return romberg(lambda u: norm(self.__call__(u, 1)), self.U[0], self.U[-1])

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
        sp = self.__call__(self.U[0])
        ep = self.__call__(self.U[-1])
        self.isClosed = (sp == ep).all()
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
        nPw = np.empty(self.Pw.shape, float)

        for i in range(self.n + 1):
            nP[i] += dv
            nPw[i] = to_homogeneous(nP[i], self.weight[i])

        self.reset(nU, nPw)

    def rotate(self, ref, ax, ang):
        """
        将曲线绕过指定点的转轴旋转一定角度
        :param ref: 参考点
        :param ax: 旋转轴方向向量
        :param ang: 旋转角
        :return: None
        """

        pass

    def to_iges(self, planar=0, periodic=0, norm_vector=np.zeros(3), form=0):
        return IGES_Entity126(self.p, self.n, planar, (1 if self.isClosed else 0), (1 if self.isPoly else 0), periodic,
                              self.U, self.weight, self.cpt, self.U[0], self.U[-1], norm_vector, form)

    def insert_knot(self, u, r=1):
        """
        插入一个节点若干次。
        :param u: 待插入节点
        :param r: 插入的次数，要求s+r<=p, 其中s为u在原节点矢量中的重复度,p为曲线次数
        :return: None
        """

        if r <= 0:
            raise ValueError('Invalid times!')

        '''Insert'''
        s = sum(x == u for x in self.U)  # Counts of duplicates
        if s + r > self.p:
            raise ValueError('Too many Knot: {}\nExisting: {}, Targeting: {} Max: {}'.format(u, s, s + r, self.p))

        k = find_span(self.n, self.p, u, self.U)
        nU = np.insert(self.U, k + 1, np.full(r, u, float))  # New knot vector
        nPw = np.zeros((self.n + r + 1, 4))  # New homogeneous control points

        '''Calculate new control points'''
        Rw = np.zeros((self.p + 1, 4))  # Holding temporary points

        '''Store unchanged control points'''
        for i in range(k - self.p + 1):
            nPw[i] = np.copy(self.Pw[i])
        for i in range(k - s, self.n + 1):
            nPw[i + r] = np.copy(self.Pw[i])
        for i in range(self.p - s + 1):
            Rw[i] = np.copy(self.Pw[k - self.p + i])

        '''Insert target knot r times'''
        L = 0
        for j in range(1, r + 1):
            L = k - self.p + j
            for i in range(self.p - j - s + 1):
                alpha = (u - self.U[L + i]) / (self.U[i + k + 1] - self.U[L + i])
                Rw[i] = alpha * Rw[i + 1] + (1.0 - alpha) * Rw[i]
            nPw[L] = np.copy(Rw[0])
            nPw[k + r - j - s] = np.copy(Rw[self.p - j - s])

        '''Load remaining control points'''
        for i in range(L + 1, k - s):
            nPw[i] = np.copy(Rw[i - L])

        '''Update'''
        self.reset(nU, nPw)

    def refine(self, X):
        """
        节点细化，插入额外的节点序列
        :param X: 待插入节点序列(已按升序排好)
        """

        if len(X) == 0:
            return

        r = len(X) - 1
        nU = np.zeros(self.m + r + 2, float)  # New knot vector
        nPw = np.zeros((self.n + r + 2, 4), float)  # New homogeneous control points

        '''Knot span'''
        a = find_span(self.n, self.p, X[0], self.U)
        b = find_span(self.n, self.p, X[r], self.U) + 1

        '''Copy unchanged control points and knots'''
        for j in range(a - self.p + 1):
            nPw[j] = np.copy(self.Pw[j])
        for j in range(b - 1, self.n + 1):
            nPw[j + r + 1] = np.copy(self.Pw[j])

        for j in range(a + 1):
            nU[j] = self.U[j]
        for j in range(b + self.p, self.m + 1):
            nU[j + r + 1] = self.U[j]

        '''Insert'''
        i = b + self.p - 1
        k = b + self.p + r
        for j in range(r, -1, -1):
            while X[j] <= self.U[i] and i > a:
                nPw[k - self.p - 1] = np.copy(self.Pw[i - self.p - 1])
                nU[k] = self.U[i]
                k -= 1
                i -= 1

            nPw[k - self.p - 1] = np.copy(nPw[k - self.p])

            for l in range(1, self.p + 1):
                index = k - self.p + l
                alpha = nU[k + l] - X[j]
                if equal(alpha, 0.0):
                    nPw[index - 1] = np.copy(nPw[index])
                else:
                    alpha /= (nU[k + l] - self.U[i - self.p + l])
                    nPw[index - 1] = alpha * nPw[index - 1] + (1.0 - alpha) * nPw[index]

            nU[k] = X[j]
            k -= 1

        self.reset(nU, nPw)

    def decompose(self, self_update=False, return_raw=False):
        """
        将NURBS曲线分解为Bezier曲线段
        :return: 分解后的节点矢量与控制点
        """

        '''New knot vector and control points'''
        k = 0
        val = np.unique(self.U)
        nU = np.empty((self.p + 1) * len(val), float)
        Qw = np.empty((len(val) - 1, self.p + 1, 4), float)

        for v in val:
            for j in range(self.p + 1):
                nU[k] = v
                k += 1

        '''Calculate new control points'''
        alphas = np.empty(self.p, float)
        a = self.p
        b = self.p + 1
        nb = 0
        for i in range(self.p + 1):
            Qw[nb][i] = np.copy(self.Pw[i])

        while b < self.m:
            i = b
            while b < self.m and equal(self.U[b + 1], self.U[b]):
                b += 1
            mult = b - i + 1
            if mult < self.p:
                numer = self.U[b] - self.U[a]
                for j in range(self.p, mult, -1):
                    alphas[j - mult - 1] = numer / (self.U[a + j] - self.U[a])
                r = self.p - mult
                for j in range(1, r + 1):
                    save = r - j
                    s = mult + j
                    for k in range(self.p, s - 1, -1):
                        alpha = alphas[k - s]
                        Qw[nb][k] = alpha * Qw[nb][k] + (1.0 - alpha) * Qw[nb][k - 1]
                    if b < self.m:
                        Qw[nb + 1][save] = np.copy(Qw[nb][self.p])

            nb += 1
            if b < self.m:
                for i in range(self.p - mult, self.p + 1):
                    Qw[nb][i] = np.copy(self.Pw[b - self.p + i])
                a = b
                b += 1

        nPw = np.empty((len(Qw) * (self.p + 1), 4), float)
        k = 0
        for i in range(len(Qw)):
            for j in range(self.p + 1):
                nPw[k] = np.copy(Qw[i][j])
                k += 1

        if self_update:
            self.reset(nU, nPw)

        if return_raw:
            return nU, nPw
        else:
            return NURBS_Curve(nU, nPw)

    def elevate(self, t: int):
        """
        将曲线升阶t次
        :param t: 升阶次数
        :return: None.
        """

        if t <= 0:
            return

        '''Basic counting'''
        val, cnt = np.unique(self.U, return_counts=True)
        for i in range(len(cnt)):
            cnt[i] += t

        mm = sum(cnt) - 1
        pp = self.p + t
        nn = mm - pp - 1
        pp2 = int(pp / 2)

        '''New knot vector and control points'''
        nU = np.zeros(mm + 1)
        Qw = np.zeros((nn + 1, 4))

        '''Local arrays'''
        bezalfs = np.zeros((pp + 1, self.p + 1))
        bpts = np.zeros((self.p + 1, 4))
        ebpts = np.zeros((pp + 1, 4))
        Nextbpts = np.zeros((self.p - 1, 4))
        alfs = np.zeros(self.p - 1)

        # 计算升阶的次数
        bezalfs[0][0] = bezalfs[pp][self.p] = 1.0

        for i in range(1, pp2 + 1):
            inv = 1.0 / comb(pp, i)
            mpi = min(self.p, i)
            for j in range(max(0, i - t), mpi + 1):
                bezalfs[i][j] = inv * comb(self.p, j) * comb(t, i - j)

        for i in range(pp2 + 1, pp):
            mpi = min(self.p, i)
            for j in range(max(0, i - t), mpi + 1):
                bezalfs[i][j] = bezalfs[pp - i][self.p - j]

        mh = pp
        kind = pp + 1
        r = -1
        a = self.p
        b = self.p + 1
        cind = 1
        ua = self.U[0]
        Qw[0] = np.copy(self.Pw[0])
        for i in range(pp + 1):
            nU[i] = ua

        # 初始化第一个bezier段
        for i in range(self.p + 1):
            bpts[i] = np.copy(self.Pw[i])

        # 沿节点矢量进行循环
        while b < self.m:
            i = b
            while b < self.m and equal(self.U[b], self.U[b + 1]):
                b += 1
            mul = b - i + 1
            mh += (mul + t)
            ub = self.U[b]
            oldr = r
            r = self.p - mul

            # 插入节点u(b) r次
            lbz = int((oldr + 2) / 2) if oldr > 0 else 1
            rbz = int(pp - (r + 1) / 2) if r > 0 else pp

            if r > 0:
                # 插入节点以获得bezier曲线段
                numer = ub - ua
                for k in range(self.p, mul, -1):
                    alfs[k - mul - 1] = numer / (self.U[a + k] - ua)
                for j in range(1, r + 1):
                    save = r - j
                    s = mul + j
                    for k in range(self.p, s - 1, -1):
                        bpts[k] = alfs[k - s] * bpts[k] + (1.0 - alfs[k - s]) * bpts[k - 1]
                    Nextbpts[save] = np.copy(bpts[self.p])

            # 对bezier曲线段升阶
            for i in range(lbz, pp + 1):
                ebpts[i].fill(0.0)
                mpi = min(self.p, i)
                for j in range(max(0, i - t), mpi + 1):
                    ebpts[i] += bezalfs[i][j] * bpts[j]

            if oldr > 1:
                first = kind - 2
                last = kind
                den = ub - ua
                bet = (ub - nU[kind - 1]) / den
                for tr in range(1, oldr):
                    i = first
                    j = last
                    kj = j - kind + 1
                    while j - i > tr:
                        if i < cind:
                            alf = (ub - nU[i]) / (ua - nU[i])
                            if i >= 9:
                                print('debug')
                            Qw[i] = alf * Qw[i] + (1.0 - alf) * Qw[i - 1]
                        if j >= lbz:
                            if j - tr <= kind - pp + oldr:
                                gam = (ub - nU[j - tr]) / den
                                ebpts[kj] = gam * ebpts[kj] + (1.0 - gam) * ebpts[kj + 1]
                            else:
                                ebpts[kj] = bet * ebpts[kj] + (1.0 - bet) * ebpts[kj + 1]
                        i += 1
                        j -= 1
                        kj -= 1
                    first -= 1
                    last += 1

            # 消去节点u=U[a]结束
            if a != self.p:
                for i in range(pp - oldr):
                    nU[kind] = ua
                    kind += 1

            # 将控制点存入Qw
            for j in range(lbz, rbz + 1):
                if cind >= 9:
                    print('debug')
                Qw[cind] = np.copy(ebpts[j])
                cind += 1

            if b < self.m:
                # 为下一次循环做准备
                for j in range(r):
                    bpts[j] = np.copy(Nextbpts[j])
                for j in range(r, self.p + 1):
                    bpts[j] = np.copy(self.Pw[b - self.p + j])
                a = b
                b += 1
                ua = ub
            else:
                for i in range(pp + 1):
                    nU[kind + i] = ub

        '''Update'''
        self.reset(nU, Qw)


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
        cdir = normalize(np.cross(norm_vector, ep - sp))
        center = 0.5 * (sp + ep) + cdir * w

        '''Rotate and pan'''
        arc = cls(radius, np.rad2deg(theta))
        base1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
        xdir = sp - center
        base2 = np.array([xdir, np.cross(norm_vector, xdir), norm_vector], float)

        mrot = DCM(base1, base2).rot_matrix
        mrot_trans = np.transpose(mrot)
        pts = np.copy(arc.cpt)
        wg = np.copy(arc.weight)
        pw = np.empty((len(pts), 4))
        for i in range(len(pts)):
            pts[i] = pts[i] * mrot_trans
            pts[i] += center
            pw[i] = to_homogeneous(pts[i], wg[i])

        '''Reconstruct'''
        arc.reset(arc.U, pw)
        return arc
