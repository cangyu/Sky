from numpy.linalg import norm
from scipy.linalg import solve
from scipy.interpolate import BSpline, make_interp_spline
from scipy.integrate import romberg
from scipy.misc import comb
from scipy.special import factorial
from src.nurbs.utility import *
from src.nurbs.transform import DCM, Quaternion
from src.iges.iges_entity110 import IGES_Entity110
from src.iges.iges_entity112 import IGES_Entity112
from src.iges.iges_entity126 import IGES_Entity126


class ClampedNURBSCrv(object):
    def __init__(self, u, pw):
        """
        NURBS曲线
        :param u:节点矢量 
        :param pw: 带权控制点
        """

        self.U = np.copy(u)
        self.Pw = np.copy(pw)

        self.spl = BSpline(self.U, self.Pw, self.p)

    @property
    def m(self):
        """
        最后一个节点下标
        """

        return len(self.U) - 1

    @property
    def n(self):
        """
        最后一个控制点下标
        """

        return len(self.Pw) - 1

    @property
    def p(self):
        """
        曲线次数
        """

        return self.m - self.n - 1

    @property
    def weight(self):
        """
        权系数
        """

        return self.Pw[:, -1]

    @property
    def cpt(self):
        """
        求控制点
        """

        ans = np.zeros((len(self.Pw), 3))
        for i in range(len(self.Pw)):
            ans[i] = to_cartesian(self.Pw[i])
        return ans

    def __call__(self, u, d=0, return_cartesian=True):
        """
        计算曲线上的点
        :param u: 目标参数
        :param d: 求导次数
        :param return_cartesian: 返回结果是否带权
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

    def reset(self, u, pw):
        """
        根据新的节点矢量和带权控制点重新初始化曲线
        :param u: New Knot vector.
        :param pw: New control points.
        :return: None.
        """

        self.U = np.copy(u)
        self.Pw = np.copy(pw)
        self.spl = BSpline(self.U, self.Pw, self.p)

    def reverse(self):
        """
        曲线反向
        """

        npw = self.Pw[::-1, :]
        nu = np.ones(self.m + 1) - self.U[::-1]
        self.reset(nu, npw)

    def pan(self, delta):
        """
        曲线整体平移
        :param delta: 偏移矢量
        :return: None.
        """

        dv = np.zeros(3)
        array_smart_copy(delta, dv)
        npw = np.empty(self.Pw.shape, float)

        for i in range(self.n + 1):
            cv = to_cartesian(self.Pw[i])
            cv += dv
            npw[i] = to_homogeneous(cv, self.Pw[i][-1])

        self.reset(self.U, npw)

    def rotate(self, ref, ax, ang):
        """
        将曲线绕过指定点的转轴旋转一定角度
        :param ref: 参考点
        :param ax: 旋转轴方向向量，按右手定则确定角度的正方向
        :param ang: 旋转角(Degree)
        :return: None
        """

        q = Quaternion.from_u_theta(ax, math.radians(ang))
        npw = np.empty(self.Pw.shape, float)
        for i in range(self.n + 1):
            cv = to_cartesian(self.Pw[i])
            cv -= ref
            cv = ref + q.rotate(cv)
            npw[i] = to_homogeneous(cv, self.Pw[i][-1])

        self.reset(self.U, npw)

    @property
    def start(self):
        """
        曲线起点
        """

        return to_cartesian(self.Pw[0])

    @property
    def end(self):
        """
        曲线终点
        """

        return to_cartesian(self.Pw[-1])

    def to_iges(self, *args, **kwargs):
        """
        将曲线以IGES标准中的第126号实体呈现
        :param args: 依次指示:isPlaner, isPeriodic, norm, form
        :param kwargs: Not used currently.
        :return: IGES_Entity126 Object.
        """

        planar = 0
        periodic = 0
        norm_vector = np.zeros(3)
        form = kwargs['form'] if 'form' in kwargs else 0
        if len(args) != 0:
            planar = args[0]
            periodic = args[1]
            norm_vector = args[2]

        w = self.weight
        cpt = self.cpt
        poly = 0 if (w != np.ones(w.shape)).any() else 1
        closed = 1 if equal(norm(self.end - self.start), 0.0) else 0

        return IGES_Entity126(self.p, self.n, planar, closed, poly, periodic, self.U, w, cpt, self.U[0], self.U[-1], norm_vector, form)

    def insert_knot(self, u, r=1):
        """
        插入一个节点若干次
        :param u: 待插入节点
        :param r: 插入的次数，要求s+r<=p, 其中s为u在原节点矢量中的重复度,p为曲线次数
        :return: None.
        """

        if r < 0:
            raise AssertionError('Invalid times!')

        if r == 0:
            return

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

    def refine(self, extra_knots):
        """
        节点细化，插入额外的节点序列
        :param extra_knots: 待插入节点序列(已按升序排好)
        """

        if len(extra_knots) == 0:
            return

        r = len(extra_knots) - 1
        nu = np.zeros(self.m + r + 2, float)  # New knot vector
        npw = np.zeros((self.n + r + 2, 4), float)  # New homogeneous control points

        '''Knot span'''
        a = find_span(self.n, self.p, extra_knots[0], self.U)
        b = find_span(self.n, self.p, extra_knots[r], self.U) + 1

        '''Copy unchanged control points and knots'''
        for j in range(a - self.p + 1):
            npw[j] = np.copy(self.Pw[j])
        for j in range(b - 1, self.n + 1):
            npw[j + r + 1] = np.copy(self.Pw[j])

        for j in range(a + 1):
            nu[j] = self.U[j]
        for j in range(b + self.p, self.m + 1):
            nu[j + r + 1] = self.U[j]

        '''Insert'''
        i = b + self.p - 1
        k = b + self.p + r
        for j in range(r, -1, -1):
            while extra_knots[j] <= self.U[i] and i > a:
                npw[k - self.p - 1] = np.copy(self.Pw[i - self.p - 1])
                nu[k] = self.U[i]
                k -= 1
                i -= 1

            npw[k - self.p - 1] = np.copy(npw[k - self.p])

            for l in range(1, self.p + 1):
                index = k - self.p + l
                alpha = nu[k + l] - extra_knots[j]
                if equal(alpha, 0.0):
                    npw[index - 1] = np.copy(npw[index])
                else:
                    alpha /= (nu[k + l] - self.U[i - self.p + l])
                    npw[index - 1] = alpha * npw[index - 1] + (1.0 - alpha) * npw[index]

            nu[k] = extra_knots[j]
            k -= 1

        self.reset(nu, npw)

    def decompose(self, self_update=False, return_raw=False):
        """
        将NURBS曲线分解为Bezier曲线段
        :return: 分解后的节点矢量与控制点
        """

        '''New knot vector and control points'''
        val = np.unique(self.U)
        nU = np.empty((self.p + 1) * len(val), float)
        Qw = np.empty((len(val) - 1, self.p + 1, 4), float)

        k = 0
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
            return ClampedNURBSCrv(nU, nPw)

    def elevate(self, t, self_update=False, return_raw=False):
        """
        将曲线升阶t次
        :param t: 升阶次数
        :type t: int
        :param self_update: 是否更新自己
        :param return_raw: 是否返回原始形式的数据(节点矢量与带权控制点)
        :return 根据不同的选项返回不同形式的结果
        """

        if t <= 0:
            return

        '''Decompose'''
        nU, nPw = self.decompose(return_raw=True)
        pp1 = self.p + 1
        k = 0
        bezier_seg = []
        while k + pp1 < len(nU):
            cbs = BezierCrv(nU[k], nU[k + pp1], self.p, nPw[k:k + pp1])
            bezier_seg.append(cbs)
            k += pp1

        '''Elevate each bezier segment'''
        bseg_num = len(bezier_seg)
        for i in range(bseg_num):
            bezier_seg[i].elevate(t)

        '''Eliminate unnecessary knots'''
        ppt = self.p + t

        nU = np.empty((bseg_num + 1) * ppt + 2, float)
        nU[0] = bezier_seg[0].a
        k = 1
        for bsg in bezier_seg:
            tmp = bsg.a
            for i in range(ppt):
                nU[k] = tmp
                k += 1

        tmp = bezier_seg[-1].b
        for i in range(ppt + 1):
            nU[k] = tmp
            k += 1

        nPw = np.empty((bseg_num * ppt + 1, 4), float)
        k = 0
        for bsg in bezier_seg:
            for i in range(bsg.n):
                nPw[k] = np.copy(bsg.Pw[i])
                k += 1

        nPw[-1] = np.copy(bezier_seg[-1].Pw[-1])

        if self_update:
            self.reset(nU, nPw)

        if return_raw:
            return nU, nPw
        else:
            return ClampedNURBSCrv(nU, nPw)

    def segment(self, a, b):
        """
        提取其中一段
        :param a: 起始节点
        :type a: float
        :param b: 终止节点
        :type b: float
        :return: 曲线上的一段
        :rtype: ClampedNURBSCrv
        """

        pass


class BezierCrv(ClampedNURBSCrv):
    def __init__(self, a, b, p, pw):
        kv = []
        for i in range(p + 1):
            kv.append(a)
        for i in range(p + 1):
            kv.append(b)

        super(BezierCrv, self).__init__(kv, pw)

    @property
    def a(self):
        return self.U[0]

    @property
    def b(self):
        return self.U[-1]

    def elevate(self, t: int):
        """
        将Bezier曲线升阶t次
        :param t: 升阶次数
        :return: None.
        """

        if t <= 0:
            return

        nh = self.p + t + 1
        npw = np.zeros((nh, 4))
        for i in range(nh):
            for j in range(max(0, i - t), min(self.p, i) + 1):
                cc = comb(self.p, j, exact=True) * comb(t, i - j, exact=True) / comb(self.p + t, i, exact=True)
                npw[i] += cc * self.Pw[j]

        kv = []
        ph = self.p + t + 1
        for i in range(ph):
            kv.append(self.a)
        for i in range(ph):
            kv.append(self.b)

        self.reset(kv, npw)


class GlobalInterpolatedCrv(ClampedNURBSCrv):
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
        kv = calc_knot_vector(param, p)
        cpt = calc_ctrl_pts(kv, p, pts, param)

        pw = np.zeros((n + 1, dim + 1))
        for i in range(0, n + 1):
            pw[i] = to_homogeneous(cpt[i])

        super(GlobalInterpolatedCrv, self).__init__(kv, pw)


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


class Line(ClampedNURBSCrv):
    def __init__(self, a, b):
        """
        两点间直线段
        :param a: 起始点坐标
        :param b: 终点坐标
        """

        u = np.array([0, 0, 1, 1])
        pw = np.array([[0, 0, 0, 1], [0, 0, 0, 1]], float)
        array_smart_copy(a, pw[0])
        array_smart_copy(b, pw[1])

        super(Line, self).__init__(u, pw)

    def length(self):
        return pnt_dist(self.start, self.end)

    def curvature(self, u):
        return 0.0

    def to_iges(self):
        return IGES_Entity110(to_cartesian(self.Pw[0]), to_cartesian(self.Pw[-1]))


class Arc(ClampedNURBSCrv):
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
        P[0] = np.array([r, 0, 0], float)
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
        pw = np.zeros((len(pts), 4))
        for i in range(len(pts)):
            pts[i] = pts[i] * mrot_trans
            pts[i] += center
            pw[i] = to_homogeneous(pts[i], wg[i])

        '''Reconstruct'''
        arc.reset(arc.U, pw)
        return arc


class Spline(object):
    def __init__(self, pts, p=3, bc_x=([(2, 0)], [(2, 0)]), bc_y=([(2, 0)], [(2, 0)]), bc_z=([(2, 0)], [(2, 0)])):
        """
        3次样条曲线,利用SciPy中的BSpline插值，
        与NURBS_Curve相比不同之处在于其节点不是Clamped,
        与GlobalInterpolatedCrv相比不同之处在于能够控制端点的高阶导数.
        """

        '''Natural coordinates'''
        self.n = len(pts)
        self.U = np.zeros(self.n, float)
        for i in range(1, self.n):
            self.U[i] = pnt_dist(pts[i], pts[i - 1]) + self.U[i - 1]
        self.n -= 1

        '''Interpolation Function'''
        self.p = int(p)
        fx = make_interp_spline(self.U, pts[:, 0], k=self.p, bc_type=bc_x)
        fy = make_interp_spline(self.U, pts[:, 1], k=self.p, bc_type=bc_y)
        fz = make_interp_spline(self.U, pts[:, 2], k=self.p, bc_type=bc_z)
        self.f = [fx, fy, fz]

    def x(self, u, d=0):
        """
        Normalized representation of x dimension
        """

        return self.f[0](u * self.U[-1], d)

    def y(self, u, d=0):
        """
        Normalized representation of y dimension
        """

        return self.f[1](u * self.U[-1], d)

    def z(self, u, d=0):
        """
        Normalized representation of z dimension
        """

        return self.f[2](u * self.U[-1], d)

    def __call__(self, u, d=0):
        """
        求指定位置上的导矢量
        :param u: 目标参数
        :param d: 求导次数
        :return: 曲线在u处的d阶导矢
        """

        return np.array([self.x(u, d), self.y(u, d), self.z(u, d)])

    def to_iges(self):
        """
        将3次样条曲线以IGES标准中第112号实体呈现
        :return: 曲线的IGES实体表示
        :rtype IGES_Entity112
        """

        if self.p != 3:
            raise AttributeError('Should be cubic curve!')

        coefficient_matrix = np.zeros((self.n + 1, 3, self.p + 1))
        for k in range(self.n + 1):
            for i in range(3):
                for j in range(self.p + 1):
                    coefficient_matrix[k][i][j] = self.f[i](self.U[k], j) / factorial(j)

        return IGES_Entity112(self.U, coefficient_matrix)
