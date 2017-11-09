import unittest
import math
from copy import deepcopy
import numpy as np
from numpy.linalg import norm
from scipy.interpolate import BSpline
from iges import Model, Entity128
from transform import Quaternion
from misc import array_smart_copy, normalize, read_airfoil_pts
from nurbs_basis import to_cartesian, to_homogeneous, point_to_line, line_intersection
from nurbs_crv import Crv, calc_pnt_param, calc_knot_vector, calc_ctrl_pts, Arc, Line, GlobalInterpolatedCrv
from wing import WingProfile

"""
Implementation of the NURBS surface.

Note:
All the NURBS notations are in the 'Clamped' format by default.
"""


class Surf(object):
    def __init__(self, u, v, pw):
        """
        Base class for NURBS Surface.
        :param u: u方向节点矢量, n+1个元素
        :param v: v方向节点矢量，m+1个元素
        :param pw: 齐次坐标序列，(n+1)x(m+1)个元素
        """

        self.U = np.copy(u)
        self.V = np.copy(v)
        self.Pw = np.copy(pw)

        self.spl = []
        for i in range(self.n + 1):
            self.spl.append(BSpline(self.V, self.Pw[i], self.q))

    @property
    def n(self):
        """
        U方向最后一个控制点下标
        """

        return self.Pw.shape[0] - 1

    @property
    def m(self):
        """
        V方向最后一个控制点下标
        """

        return self.Pw.shape[1] - 1

    @property
    def p(self):
        """
        :return: Degree in U direction.
        :rtype: int
        """

        return len(self.U) - self.n - 2

    @property
    def q(self):
        """
        :return: Degree in V direction.
        :rtype: int
        """

        return len(self.V) - self.m - 2

    @property
    def weight(self):
        """
        权系数
        """

        return self.Pw[:, :, -1]

    @property
    def cpt(self):
        """
        Control points.
        """

        ans = np.zeros((self.n + 1, self.m + 1, 3))
        for i in range(self.n + 1):
            for j in range(self.m + 1):
                ans[i][j] = to_cartesian(self.Pw[i][j])
        return ans

    def __repr__(self):
        return 'U Knot:\n{}\nV Knot:\n{}\nControl points:\n{}'.format(self.U, self.V, self.Pw)

    def __str__(self):
        ret = 'Clamped NURBS Surface\nDegree:({},{})\n'.format(self.p, self.q)
        ret += self.__repr__()
        return ret

    def __call__(self, u, v, k=0, l=0):
        """
        求在给定位置(u,v)处的导矢量
        :param u: U方向参数
        :param v: V方向参数
        :param k: U方向求导次数
        :param l: V方向求导次数
        :return: (u,v)处偏导矢量
        """

        r = []
        for spl in self.spl:
            r.append(spl(v, l))

        rw = np.copy(r)
        spl = BSpline.construct_fast(self.U, rw, self.p)
        pw = spl(u, k)
        return to_cartesian(pw)

    def reset(self, u, v, pw):
        """
        重置曲面
        :param u: u方向节点矢量, n+1个元素
        :param v: v方向节点矢量，m+1个元素
        :param pw: 齐次坐标序列，(n+1)x(m+1)个元素
        """

        self.U = np.copy(u)
        self.V = np.copy(v)
        self.Pw = np.copy(pw)

        self.spl = []
        q = self.q
        for i in range(self.n + 1):
            self.spl.append(BSpline(self.V, self.Pw[i], q))

    def reverse(self, direction):
        """
        曲面反向
        """

        if direction not in ('U', 'V', 'UV'):
            raise ValueError('Invalid direction choice!')

        if direction in ('U', 'UV'):
            self.U = np.full(self.U.shape, self.U[0] + self.U[-1]) - self.U[::-1]
            self.Pw = self.Pw[::-1, :, :]

        if direction in ('V', 'UV'):
            self.V = np.full(self.V.shape, self.V[0] + self.V[-1]) - self.V[::-1]
            self.Pw = self.Pw[:, ::-1, :]

        self.reset(self.U, self.V, self.Pw)

    def swap(self):
        """
        交换UV方向节点矢量与控制点
        :return: None.
        """

        tmp = self.U[:]
        self.U = self.V[:]
        self.V = tmp
        self.Pw = np.transpose(self.Pw, (1, 0, 2))
        self.reset(self.U, self.V, self.Pw)

    def pan(self, delta):
        """
        曲面整体平移
        :param delta: 偏移矢量
        :return: None.
        """

        dv = np.zeros(3)
        array_smart_copy(delta, dv)

        for i in range(self.n + 1):
            for j in range(self.m + 1):
                cv = to_cartesian(self.Pw[i][j]) + dv
                self.Pw[i][j] = to_homogeneous(cv, self.Pw[i][j][-1])

        self.reset(self.U, self.V, self.Pw)

    def rotate(self, ref, ax, ang):
        """
        将曲面绕过指定点的转轴旋转一定角度
        :param ref: 参考点
        :param ax: 旋转轴方向向量，按右手定则确定角度的正方向
        :param ang: 旋转角(Degree)
        :return: None.
        """

        q = Quaternion.from_u_theta(ax, math.radians(ang))
        for i in range(self.n + 1):
            for j in range(self.m + 1):
                cv = to_cartesian(self.Pw[i][j]) - ref
                cv = ref + q.rotate(cv)
                self.Pw[i][j] = to_homogeneous(cv, self.Pw[i][j][-1])
        self.reset(self.U, self.V, self.Pw)

    def mirror(self, axis):
        """
        Mirror the surface along specified axis.
        :param axis: Direction axis.
        :return: None.
        """

        '''Defensive check'''
        if axis in ('X', 'x'):
            idx = 0
        elif axis in ('Y', 'y'):
            idx = 1
        elif axis in ('Z', 'z'):
            idx = 2
        else:
            raise ValueError("Invalid axis")

        '''Modify control points'''
        for i in range(self.n + 1):
            for j in range(self.m + 1):
                self.Pw[i][j][idx] *= -1

        '''Update'''
        self.reset(self.U, self.V, self.Pw)

    def to_iges(self, *args, **kwargs):
        """
        将曲面以IGES标准中第128号实体呈现
        :return: IGES_Entity128 Object
        """

        w = self.weight
        poly = 0 if (w != np.ones(w.shape)).any() else 1
        cpt = self.cpt
        closed_u = 0
        closed_v = 0
        periodic_u = 0
        periodic_v = 0
        form = 0

        return Entity128(self.U, self.V, self.p, self.q, self.n, self.m, cpt, w, closed_u, closed_v, poly, periodic_u, periodic_v, self.U[0], self.U[-1], self.V[0], self.V[-1], form)

    def insert_knot(self, uv, r=1, direction='U'):
        """
        曲面插入节点
        :param uv: 待插入节点值
        :param r: 插入次数
        :param direction: 插入的方向
        :return: None
        """

        if direction not in ('U', 'V'):
            raise AssertionError('Invalid direction choice!')

        if direction == 'U':
            crv_list = []
            npw = np.zeros((self.n + 2, self.m + 1, 4))
            for j in range(self.m + 1):
                cc = Crv(self.U, self.Pw[:, j])
                cc.insert_knot(uv, r)
                crv_list.append(cc)
                for i in range(self.n + 2):
                    npw[i][j] = np.copy(cc.Pw[i])

            self.reset(crv_list[0].U, self.V, npw)

        else:
            crv_list = []
            npw = np.zeros((self.n + 1, self.m + 2, 4))
            for i in range(self.n + 1):
                cc = Crv(self.V, self.Pw[i, :])
                cc.insert_knot(uv, r)
                crv_list.append(cc)
                for j in range(self.m + 2):
                    npw[i][j] = np.copy(cc.Pw[j])

            self.reset(self.U, crv_list[0].U, npw)

    def extract(self, direction, uv):
        """
        提取等参数线
        :param direction: 方向
        :param uv: 等参数值
        :return: 给定方向上的等参数线
        :rtype: Crv
        """

        if direction not in ('U', 'V'):
            raise AssertionError('Invalid direction choice!')
        if np.less(uv, 0) or np.greater(uv, 1):
            raise AssertionError('Invalid parameter!')

        if direction == 'U':
            nqw = np.zeros((self.m + 1, 4))
            for j in range(self.m + 1):
                spl = BSpline(self.U, self.Pw[:, j, :], self.p)
                nqw[j] = spl(uv)

            return Crv(self.V, nqw)

        else:
            npw = np.zeros((self.n + 1, 4))
            for i in range(self.n + 1):
                spl = BSpline(self.V, self.Pw[i, :, :], self.q)
                npw[i] = spl(uv)

            return Crv(self.U, npw)

    def refine(self, direction, extra_knot):
        """
        细化节点矢量
        :param direction: 方向选择
        :param extra_knot: 待插入节点数组
        :return: None
        """

        if direction not in ('U', 'V'):
            raise AssertionError('Invalid direction choice!')
        if len(extra_knot) == 0:
            return

        crv_list = []
        if direction == 'U':
            nh = self.n + 1 + len(extra_knot)
            npw = np.zeros((nh, self.m + 1, 4))
            for j in range(self.m + 1):
                cc = Crv(self.U, self.Pw[:, j, :])
                cc.refine(extra_knot)
                crv_list.append(cc)
                for i in range(nh):
                    npw[i][j] = np.copy(cc.Pw[i])

            self.reset(crv_list[0].U, self.V, npw)

        else:
            mh = self.m + 1 + len(extra_knot)
            npw = np.zeros((self.n + 1, mh, 4))
            for i in range(self.n + 1):
                cc = Crv(self.V, self.Pw[i, :, :])
                cc.refine(extra_knot)
                crv_list.append(cc)
                for j in range(mh):
                    npw[i][j] = np.copy(cc.Pw[j])

            self.reset(self.U, crv_list[0].U, npw)

    def elevate(self, tu=0, tv=0):
        """
        曲面升阶
        :param tu: U方向升阶次数
        :type tu: int
        :param tv: V方向升阶次数
        :type tv: int
        :return: None.
        """

        if tu < 0 or tv < 0:
            raise AssertionError('Invalid promotion!')

        if tu > 0:
            crv_list = []
            for j in range(self.m + 1):
                cc = Crv(self.U, self.Pw[:, j, :])
                cc.elevate(tu)
                crv_list.append(cc)

            nh = len(crv_list[0].Pw)
            npw = np.zeros((nh, self.m + 1, 4))
            for j in range(self.m + 1):
                for i in range(nh):
                    npw[i][j] = np.copy(crv_list[j].Pw[i])

            self.reset(crv_list[0].U, self.V, npw)

        if tv > 0:
            crv_list = []
            for i in range(self.n + 1):
                cc = Crv(self.V, self.Pw[i, :, :])
                cc.elevate(tv)
                crv_list.append(cc)

            mh = len(crv_list[0].Pw)
            npw = np.zeros((self.n + 1, mh, 4))
            for i in range(self.n + 1):
                for j in range(mh):
                    npw[i][j] = np.copy(crv_list[i].Pw[j])

            self.reset(self.U, crv_list[0].U, npw)

    def reparameterization(self, alpha, beta, gamma, delta, direction):
        """
        使用齐次有理函数将曲面重新参数化
             alpha * u + beta
        s = -------------------
             gamma * u + delta
        :param alpha: parameter
        :type alpha: float
        :param beta: parameter
        :type beta: float
        :param gamma: parameter
        :type gamma: float
        :param delta: parameter
        :type delta: float
        :param direction: parameter
        :type direction: str
        :return: None
        """

        if direction not in ('U', 'u', 'V', 'v'):
            raise AssertionError("Invalid direction parameter.")

        def g(x):
            return (alpha * x + beta) / (gamma * x + delta)

        def nbla(x):
            return gamma * x - alpha

        cpt = self.cpt
        npw = np.empty_like(self.Pw)
        wb = self.weight
        factor = 1.0

        if direction in ('U', 'u'):
            s = np.copy(list(map(g, self.U)))
            for k in range(self.p):
                factor *= nbla(s[k])
            for i in range(self.n + 1):
                factor /= nbla(s[i])
                factor *= nbla(s[i + self.p])
                for j in range(self.m + 1):
                    wb[i][j] *= factor
                    npw[i][j] = to_homogeneous(cpt[i][j], wb[i][j])
            self.reset(s, self.V, npw)
        else:
            t = np.copy(list(map(g, self.V)))
            for k in range(self.q):
                factor *= nbla(t[k])
            for j in range(self.m + 1):
                factor /= nbla(t[j])
                factor *= nbla(t[j + self.q])
                for i in range(self.n + 1):
                    wb[i][j] *= factor
                    npw[i][j] = to_homogeneous(cpt[i][j], wb[i][j])
            self.reset(self.U, t, npw)

    def standard_reparameterization(self):
        """
        将U,V两个方向的节点都统一到[0, 1]
        :return: None
        """

        '''U direction reparameterization'''
        a = self.U[0]
        b = self.U[-1]
        alpha = 1.0
        beta = -a
        gamma = 0.0
        delta = b - a
        self.reparameterization(alpha, beta, gamma, delta, 'U')

        '''V direction reparameterization'''
        a = self.V[0]
        b = self.V[-1]
        alpha = 1.0
        beta = -a
        gamma = 0.0
        delta = b - a
        self.reparameterization(alpha, beta, gamma, delta, 'V')

    @classmethod
    def split(cls, surf, u_brk, v_brk):
        """
        将曲面分割成若干子部分
        :param surf: Surface to be split
        :type surf: Surf
        :param u_brk: breaking knot in u-direction
        :param v_brk: breaking knot in v-direction
        :return: Collection of split surf
        """

        cp = surf.p
        cq = surf.q

        '''Pre-check'''
        if len(u_brk) != 0 and (min(u_brk) <= 0 or max(u_brk) >= 1):
            raise AssertionError("Invalid input.")

        if len(v_brk) != 0 and (min(v_brk) <= 0 or max(v_brk) >= 1):
            raise AssertionError("Invalid input.")

        '''Copy back break knot info'''
        uspk = sorted(u_brk)
        uspk.append(1.0)
        vspk = sorted(v_brk)
        vspk.append(1.0)

        '''Statistic current surf knot info'''
        uval, ucnt = np.unique(surf.U, return_counts=True)
        ukdt = dict(zip(uval, ucnt))

        vval, vcnt = np.unique(surf.V, return_counts=True)
        vkdt = dict(zip(vval, vcnt))

        '''Construct knot to be inserted'''
        uek = []
        for u in uspk:
            exist_cnt = ukdt.get(u) if u in ukdt else 0
            tc = cp - exist_cnt
            if tc > 0:
                for k in range(tc):
                    uek.append(u)

        vek = []
        for v in vspk:
            exist_cnt = vkdt.get(v) if v in vkdt else 0
            tc = cq - exist_cnt
            if tc > 0:
                for k in range(tc):
                    vek.append(v)

        '''Insert knots'''
        vsrf = deepcopy(surf)
        vsrf.refine('U', np.copy(uek))
        vsrf.refine('V', np.copy(vek))

        '''Build knot segment'''
        usdt = []
        uprev = 0.0
        ucki = cp + 1
        for u in uspk:
            cu = []
            for k in range(cp + 1):
                cu.append(uprev)
            while ucki < len(vsrf.U) and vsrf.U[ucki] <= u:
                cu.append(vsrf.U[ucki])
                ucki += 1
            if ucki < len(vsrf.U):
                cu.append(u)
            uprev = u
            usdt.append(cu)

        vsdt = []
        vprev = 0.0
        vcki = cq + 1
        for v in vspk:
            cv = []
            for k in range(cq + 1):
                cv.append(vprev)
            while vcki < len(vsrf.V) and vsrf.V[vcki] <= v:
                cv.append(vsrf.V[vcki])
                vcki += 1
            if vcki < len(vsrf.V):
                cv.append(v)
            vprev = v
            vsdt.append(cv)

        '''Extract control points'''
        ret = []
        ucpis = 0
        for useg in usdt:
            ucpn = len(useg) - cp - 1
            vcpis = 0
            csrf_seg = []
            for vseg in vsdt:
                vcpn = len(vseg) - cq - 1
                cpt = vsrf.Pw[ucpis:ucpis + ucpn, vcpis:vcpis + vcpn]
                vcpis += vcpn - 1
                csrf = Surf(useg, vseg, cpt)
                csrf.standard_reparameterization()
                csrf_seg.append(csrf)
            ucpis += ucpn - 1
            ret.append(csrf_seg)

        return ret


class GlobalInterpolatedSurf(Surf):
    def __init__(self, pts, p, q, u_method='centripetal', v_method='chord'):
        """
        (n+1)x(m+1)个数据点全局插值，非有理
        不能很好处理局部数据点共面，需小心使用
        :param pts: 待插值数据点
        :param p: u方向次数
        :param q: v方向次数
        :param u_method: u方向参数计算方法
        :param v_method: v方向参数计算方法
        """

        n, m, dim = pts.shape
        n -= 1
        m -= 1

        U = np.zeros(n + 1)
        V = np.zeros(m + 1)
        U[-1] = 1.0
        V[-1] = 1.0

        '''Parameters of U direction'''
        dist = np.zeros((n + 1, m + 1))
        for j in range(0, m + 1):
            td = calc_pnt_param(pts[:, j], u_method)
            for i in range(0, n + 1):
                dist[i][j] = td[i]
        for i in range(0, n):
            U[i] = np.mean(dist[i])

        '''Parameters of V Direction'''
        for i in range(0, n + 1):
            td = calc_pnt_param(pts[i], v_method)
            for j in range(0, m + 1):
                dist[i][j] = td[j]
        for j in range(0, m):
            V[j] = np.mean(dist[:, j])

        '''Knot Vectors'''
        u_knot = calc_knot_vector(U, p)
        v_knot = calc_knot_vector(V, q)

        '''Control Points'''
        R = np.zeros((n + 1, m + 1, dim))
        for j in range(0, m + 1):
            tp = calc_ctrl_pts(u_knot, p, pts[:, j], U)
            for i in range(0, n + 1):
                R[i][j] = tp[i]

        P = np.zeros((n + 1, m + 1, dim))
        for i in range(0, n + 1):
            P[i] = calc_ctrl_pts(v_knot, q, R[i], V)

        Pw = np.zeros((n + 1, m + 1, dim + 1))
        for i in range(0, n + 1):
            for j in range(0, m + 1):
                Pw[i][j] = to_homogeneous(P[i][j])

        super(GlobalInterpolatedSurf, self).__init__(u_knot, v_knot, Pw)


class BilinearSurf(Surf):
    def __init__(self, p):
        """
        双线性曲面

        ^ V direction
        |
        |P[0][1]        P[1][1]
        ----------------
        |              |
        |              |
        |     SURF     |
        |              |
        |P[0][0]       |P[1][0]
        --------------------------> U direction

        :param p:4个角点, 2x2
        """

        u_vec = np.array([0, 0, 1, 1], float)
        v_vec = np.array([0, 0, 1, 1], float)

        ul, vl, dim = p.shape
        assert ul == 2 and vl == 2

        pw = np.ones((ul, vl, 4), float)
        for i in range(ul):
            for j in range(vl):
                for d in range(dim):
                    pw[i][j][d] = p[i][j][d]

        super(BilinearSurf, self).__init__(u_vec, v_vec, pw)


class ExtrudedSurf(Surf):
    def __init__(self, crv, direction):
        """
        拉伸曲面
        :param crv: Curve to be extruded.
        :type crv: Crv
        :param direction: Direction vector.
        """

        U = np.copy(crv.U)
        V = np.array([0, 0, 1, 1], float)
        n = len(crv.cpt)
        Pw = np.zeros((n, 2, 4))
        for i in range(n):
            Pw[i][0] = Pw[i][1] = np.copy(crv.Pw[i])
            wdir = to_homogeneous(direction, Pw[i][0][3])
            for d in range(3):
                Pw[i][1][d] += wdir[d]

        super(ExtrudedSurf, self).__init__(U, V, Pw)


class RuledSurf(Surf):
    def __init__(self, _c1, _c2):
        """
        生成V方向的直纹面,即两条曲线之间的线性插值
        :param _c1: 第1条曲线
        :type _c1: Crv
        :param _c2: 第2条曲线
        :type _c2: ClampedNURBSCrv
        """

        '''Not change original curve'''
        c1 = deepcopy(_c1)
        c2 = deepcopy(_c2)

        '''Check'''
        if not math.isclose(c1.U[0], c2.U[0]):
            raise ValueError('Incompatible starting knot!')
        if not math.isclose(c1.U[-1], c2.U[-1]):
            raise ValueError('Incompatible ending knot!')

        '''Knot vector'''
        p = max(c1.p, c2.p)
        c1.elevate(p - c1.p)
        c2.elevate(p - c2.p)

        if len(c1.U) != len(c2.U) or not math.isclose(norm(c1.U - c2.U), 0):
            all_knot = merge_knot(c1.U, c2.U)
            x1 = different_knot(all_knot, c1.U)
            x2 = different_knot(all_knot, c2.U)
            c1.refine(x1)
            c2.refine(x2)

        uknot = c1.U
        vknot = np.array([0, 0, 1, 1], float)

        '''Control points'''
        pw = np.zeros((len(c1.Pw), 2, 4))
        for i in range(len(c1.Pw)):
            pw[i][0] = np.copy(c1.Pw[i])
            pw[i][1] = np.copy(c2.Pw[i])

        super(RuledSurf, self).__init__(uknot, vknot, pw)


class RevolvedSurf(Surf):
    def __init__(self, center, axis, theta, crv):
        """
        曲线绕过指定点的轴线旋转指定角度得到的曲面
        :param center: 旋转中心
        :param axis: 旋转轴，正方向按右手法则给定
        :param theta: 旋转角度
        :type theta: float
        :param crv: 母线
        :type crv: Crv
        """

        while theta <= 0:
            theta += 360
        while theta > 360:
            theta -= 360

        '''Basic variables'''
        narcs = int(math.ceil(theta / 90))
        theta = math.radians(theta)
        delta_theta = theta / narcs
        wm = math.cos(delta_theta / 2)

        '''U Direction knots'''
        u_knot = np.zeros(2 * narcs + 4)
        u_knot[-1] = u_knot[-2] = u_knot[-3] = 1.0
        delta_knot = 1.0 / narcs
        for i in range(1, narcs):
            cur_index = 1 + 2 * i
            u_knot[cur_index] = u_knot[cur_index + 1] = i * delta_knot

        '''Pre-compute sine and cosine stuff'''
        sines = np.zeros(narcs + 1, float)
        cosines = np.ones(narcs + 1, float)
        angle = 0.0
        for i in range(1, narcs + 1):
            angle += delta_theta
            cosines[i] = math.cos(angle)
            sines[i] = math.sin(angle)

        '''Compute weighted control points on each line'''
        Pj = crv.cpt
        wj = crv.weight
        m = crv.n
        npw = np.zeros((len(u_knot) - 2 - 1, m + 1, 4), float)
        for j in range(m + 1):
            O = point_to_line(Pj[j], center, axis)
            X = Pj[j] - O
            r = norm(X)
            X = normalize(X)
            Y = np.cross(axis, X)
            npw[0][j] = crv.Pw[j]
            P0 = Pj[j]
            T0 = Y
            index = 0
            for i in range(1, narcs + 1):
                P2 = O + r * cosines[i] * X + r * sines[i] * Y
                npw[index + 2][j] = to_homogeneous(P2, wj[j])
                T2 = -sines[i] * X + cosines[i] * Y
                npw[index + 1][j] = to_homogeneous(line_intersection(P0, T0, P2, T2), wm * wj[j])
                index += 2
                if i < narcs:
                    P0 = P2
                    T0 = T2

        super(RevolvedSurf, self).__init__(u_knot, crv.U, npw)


class Coons(Surf):
    def __init__(self, c0u, c1u, c0v, c1v):
        """
        双线性混合Coons曲面

         ^ V direction
         |
         |     c1u
         ------->--------
         |              |
         |              |
     c0v ^     SURF     ^ c1v
         |              |
         |              |
         ------->-----------> U direction
               c0u

        :param c0u:沿U方向第1条曲线
        :type c0u: ClampedNURBSCrv
        :param c1u:沿U方向第2条曲线
        :type c1u: ClampedNURBSCrv
        :param c0v:沿V方向第1条曲线
        :type c0v: ClampedNURBSCrv
        :param c1v:沿V方向第2条曲线
        :type c1v: Crv
        """

        '''Check 4 corners'''
        assert math.isclose(norm(c0u(0) - c0v(0)), 0.0)
        assert math.isclose(norm(c0u(1) - c1v(0)), 0.0)
        assert math.isclose(norm(c1v(1) - c1u(1)), 0.0)
        assert math.isclose(norm(c0v(1) - c1u(0)), 0.0)

        '''Corner points'''
        s = np.zeros((2, 2, 3))
        s[0][0] = np.copy(c0u(0))
        s[0][1] = np.copy(c0v(1))
        s[1][0] = np.copy(c1v(0))
        s[1][1] = np.copy(c1u(1))

        '''Base surf'''
        r1 = RuledSurf(c0u, c1u)
        r2 = RuledSurf(c0v, c1v)
        r2.swap()
        t = BilinearSurf(s)

        '''Elevate to same order'''
        pu = max(r1.p, r2.p, t.p)
        pv = max(r1.q, r2.q, t.q)
        r1.elevate(pu - r1.p, pv - r1.q)
        r2.elevate(pu - r2.p, pv - r2.q)
        t.elevate(pu - t.p, pv - t.q)

        '''Unify knot vector'''
        xu = merge_knot(merge_knot(r1.U, r2.U), t.U)
        xv = merge_knot(merge_knot(r1.V, r2.V), t.V)

        xr1u = different_knot(xu, r1.U)
        xr2u = different_knot(xu, r2.U)
        xtu = different_knot(xu, t.U)

        xr1v = different_knot(xv, r1.V)
        xr2v = different_knot(xv, r2.V)
        xtv = different_knot(xv, t.V)

        r1.refine('U', xr1u)
        r1.refine('V', xr1v)

        r2.refine('U', xr2u)
        r2.refine('V', xr2v)

        t.refine('U', xtu)
        t.refine('V', xtv)

        '''Calculate new control points'''
        npw = r1.Pw + r2.Pw - t.Pw

        super(Coons, self).__init__(xu, xv, npw)


class Skinned(Surf):
    def __init__(self, crv, p, q, v_method='chord'):
        """
        蒙皮曲面，非有理
        :param crv: 非有理曲线集合
        :param p: 目标曲面u方向次数(曲线方向)
        :param q: 目标曲面v方向次数(展向)
        :param v_method: v方向插值方法
        """

        '''Promote all curve to p order'''
        crv_list = []
        for c in crv:
            cc = deepcopy(c)
            cc.elevate(p - cc.p)
            crv_list.append(cc)

        '''Merge all knot vectors in U direction'''
        u_knot = np.copy(crv_list[0].U)
        for i in range(1, len(crv_list)):
            u_knot = merge_knot(u_knot, crv_list[i].U)

        '''Unify all curve knot vector'''
        for c in crv_list:
            xu = different_knot(u_knot, c.U)
            c.refine(xu)

        '''Knot vector in V direction'''
        n = len(u_knot) - p - 2
        m = len(crv_list) - 1
        pnt = np.zeros((n + 1, m + 1, 3))
        for j in range(m + 1):
            for i in range(n + 1):
                pnt[i][j] = to_cartesian(crv_list[j].Pw[i])

        v_param = np.zeros((n + 1, m + 1))
        vk = []
        for i in range(n + 1):
            v_param[i] = calc_pnt_param(pnt[i], v_method)
            vk.append(calc_knot_vector(v_param[i], q))

        v_knot = np.mean(vk, axis=0)

        '''Calculate control points'''
        Q = np.zeros((n + 1, m + 1, 3))
        Qw = np.zeros((n + 1, m + 1, 4))

        for i in range(n + 1):
            Q[i] = calc_ctrl_pts(v_knot, q, pnt[i], v_param[i])

        for i in range(n + 1):
            for j in range(m + 1):
                Qw[i][j] = to_homogeneous(Q[i][j])

        super(Skinned, self).__init__(u_knot, v_knot, Qw)


def merge_knot(lhs, rhs):
    """
    合并两个节点矢量
    :param lhs: 第1个节点矢量
    :param rhs: 第2个节点矢量
    :return: 合并后的节点矢量, lhs union rhs
    """

    lval, lcnt = np.unique(lhs, return_counts=True)
    rval, rcnt = np.unique(rhs, return_counts=True)
    val = np.unique(np.concatenate((lval, rval)))

    ans = []
    for v in val:
        if v in lval and v in rval:
            li = np.searchsorted(lval, v)
            ri = np.searchsorted(rval, v)
            cc = max(lcnt[li], rcnt[ri])
            for i in range(cc):
                ans.append(v)
        else:
            if v in lval:
                li = np.searchsorted(lval, v)
                for i in range(lcnt[li]):
                    ans.append(v)
            else:
                ri = np.searchsorted(rval, v)
                for i in range(rcnt[ri]):
                    ans.append(v)

    return np.copy(ans)


def different_knot(lhs, rhs):
    """
    求两个节点矢量中的不同部分
    :param lhs: 第1个节点矢量
    :param rhs: 第2个节点矢量
    :return: lhs subtract rhs
    """

    lval, lcnt = np.unique(lhs, return_counts=True)
    rval, rcnt = np.unique(rhs, return_counts=True)

    '''Count each'''
    val = []
    cnt = []
    for i in range(0, len(lval)):
        if lval[i] in rval:
            k = np.searchsorted(rval, lval[i])
            lvc = lcnt[i]
            rvc = rcnt[k]
            if lvc > rvc:
                val.append(lval[i])
                cnt.append(lvc - rvc)
        else:
            val.append(lval[i])
            cnt.append(lcnt[i])

    '''Assemble'''
    ans = np.zeros(int(sum(cnt)))
    k = 0
    for i in range(0, len(val)):
        for j in range(0, cnt[i]):
            ans[k] = val[i]
            k += 1

    return ans


class NURBSSurfTester(unittest.TestCase):
    def test_call(self):
        pass

    def test_bilinear(self):
        l = 10.0
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        n = len(p)
        iges_model = Model()
        for k in range(n):
            s = BilinearSurf(p[k])
            iges_model.add(s.to_iges())
        iges_model.save('test_bilinear.igs')
        self.assertTrue(True)

    def test_pan(self):
        l = 10.0
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])
        d = [(0, 0, 10), (0, 0, 10)]

        self.assertTrue(len(p) == len(d))
        n = len(p)

        iges_model = Model()
        for k in range(n):
            iges_model.clear()
            s = BilinearSurf(p[k])
            iges_model.add(s.to_iges())
            s.pan(d[k])
            iges_model.add(s.to_iges())
            iges_model.save('test_pan-{}.igs'.format(k))

    def test_rotate(self):
        l = 10.0
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])
        ref = [[l, l, l], [l, l, l]]
        axis = [[0, 0, 5], [0, 0, 5]]
        ang = [45, 45]

        n = len(p)
        iges_model = Model()
        for k in range(n):
            iges_model.clear()
            s = BilinearSurf(p[k])
            iges_model.add(s.to_iges())
            s.rotate(ref[k], axis[k], ang[k])
            iges_model.add(s.to_iges())
            iges_model.save('test_rotate-{}.igs'.format(k))
        self.assertTrue(True)

    def test_reverse(self):
        l = 12.34
        p = np.array([[[[0, 0, 0], [0, l, l]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]],
                       [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]],
                       [[l, 0, l], [l, l, 0]]]])
        axis = ['U', 'V', 'UV', 'U', 'V', 'UV']

        iges_model = Model()
        for k, bp in enumerate(p):
            iges_model.clear()
            s = BilinearSurf(bp)
            iges_model.add(s.to_iges())
            s.reverse(axis[k])
            iges_model.add(s.to_iges())
            iges_model.save('test_reverse-{}.igs'.format(k))

        self.assertTrue(True)

    def test_swap(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        for k, bpt in enumerate(p):
            s = BilinearSurf(bpt)
            ss = deepcopy(s)
            s.elevate(1, 2)
            ss.elevate(1, 2)
            s.swap()
            iges_model.clear()
            iges_model.add(s.to_iges())
            iges_model.add(ss.to_iges())
            iges_model.save('test_swap-{}.igs'.format(k))

            ni, nj = 50, 30
            u_dist, v_dist = np.meshgrid(np.linspace(0, 1, ni), np.linspace(0, 1, nj), indexing='ij')
            ps = np.zeros((ni, nj, 3))
            pss = np.zeros((ni, nj, 3))

            for i in range(ni):
                for j in range(nj):
                    ps[i][j] = s(v_dist[i][j], u_dist[i][j])
                    pss[i][j] = ss(u_dist[i][j], v_dist[i][j])

            self.assertTrue(np.allclose(ps, pss))

    def test_insert(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        print('Test insert utility.')
        for k, dt in enumerate(p):
            print('Test-{}'.format(k))
            iges_model.clear()
            s = BilinearSurf(dt)
            print('Before:{}'.format(s))
            iges_model.add(s.to_iges())
            s.insert_knot(0.5, 1, 'U')
            s.insert_knot(0.5, 1, 'V')
            print('After:{}'.format(s))
            iges_model.add(s.to_iges())
            iges_model.save('test_insert-{}.igs'.format(k))

        self.assertTrue(True)

    def test_refine(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        print('Test refine utility.')
        for k, dt in enumerate(p):
            print('Test-{}'.format(k))
            iges_model.clear()
            s = BilinearSurf(dt)
            ss = deepcopy(s)
            print('Before:{}'.format(s))
            iges_model.add(s.to_iges())
            s.refine('U', [0.2, 0.3, 0.6, 0.7])
            s.refine('V', [0.1, 0.4, 0.8, 0.9])
            print('After:{}'.format(s))
            iges_model.add(s.to_iges())
            iges_model.save('test_refine-{}.igs'.format(k))

            ni, nj = 50, 30
            u_dist = np.linspace(0, 1, ni)
            v_dist = np.linspace(0, 1, nj)
            ps = np.zeros((ni, nj, 3))
            pss = np.zeros((ni, nj, 3))

            for i in range(ni):
                for j in range(nj):
                    ps[i][j] = s(u_dist[i], v_dist[j])
                    pss[i][j] = ss(u_dist[i], v_dist[j])

            self.assertTrue(np.allclose(ps, pss))

    def test_elevate(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])

        iges_model = Model()
        print('Test elevate utility.')
        for k, dt in enumerate(p):
            print('Test-{}'.format(k))
            iges_model.clear()
            s = BilinearSurf(dt)
            ss = deepcopy(s)
            print('Before:{}'.format(s))
            iges_model.add(s.to_iges())
            s.elevate(1, 2)
            print('After:{}'.format(s))
            iges_model.add(s.to_iges())
            iges_model.save('test_elevate-{}.igs'.format(k))

            ni, nj = 50, 30
            u_dist = np.linspace(0, 1, ni)
            v_dist = np.linspace(0, 1, nj)
            ps = np.zeros((ni, nj, 3))
            pss = np.zeros((ni, nj, 3))
            for i in range(ni):
                for j in range(nj):
                    ps[i][j] = s(u_dist[i], v_dist[j])
                    pss[i][j] = ss(u_dist[i], v_dist[j])

            self.assertTrue(np.allclose(ps, pss))

    def test_split(self):
        print('Test split utility.')
        iges_model = Model()
        s = Surf(np.array([0, 0, 0, 0, 1., 1., 1., 1.]),
                 np.array([0, 0, 0, 0.5, 1, 1, 1]),
                 np.array([[[0, 0, 0, 1.], [0, 1, 1, 1], [0, 2, 3, 1], [0, 3, 2, 1]],
                           [[1, 0, 0, 1.], [1, 1, 2, 1], [1, 3, 5, 1], [1, 4, 2, 1]],
                           [[2, 0, 0, 1.], [2, 2, 2, 1], [2, 3, 6, 1], [2, 5, 7, 1]],
                           [[3, 0, 0, 1.], [3, 1, 1, 1], [3, 2, 2, 1], [3, 3, 3, 1]]]))
        iges_model.add(s.to_iges())
        print('Original:{}'.format(s))
        ss = Surf.split(s, (0.4, 0.5, 0.6), [0.2, 0.5, 0.7])
        k = 0
        for ue in ss:
            for ve in ue:
                iges_model.add(ve.to_iges())
                print('Segment {}:{}'.format(k, ve))
                k += 1
        iges_model.save('test_split.igs')
        self.assertTrue(True)

    def test_extract(self):
        l = 11.11
        p = np.array([[[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, 0], [0, l, l]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]],
                      [[[0, 0, l], [0, l, 0]], [[l, 0, l], [l, l, 0]]]])
        ext = [('U', 0), ('U', 0.5), ('U', 1), ('V', 0), ('V', 0.5), ('V', 1),
               ('U', 0), ('U', 0.5), ('U', 1), ('V', 0), ('V', 0.5), ('V', 1)]
        self.assertTrue(len(p) == len(ext))

        iges_model = Model()
        print('Test extract utility.')
        for k, dt in enumerate(p):
            print('Case-{}'.format(k))
            s = BilinearSurf(dt)
            print('Surf:\n{}'.format(s))
            mid = s.extract(ext[k][0], ext[k][1])
            print('Curve at {}={}:\n{}'.format(ext[k][0], ext[k][1], mid))
            iges_model.clear()
            iges_model.add(s.to_iges())
            iges_model.add(mid.to_iges())
            iges_model.save('test_extract-{}_{}={}.igs'.format(k, ext[k][0], ext[k][1]))

    def test_extrude(self):
        crv = [Arc.from_2pnt((30, 45, 69), (44, 66, 88), 75, (1, -1, 1)),
               Arc.from_2pnt((0, 50, 0), (50, 0, 0), 315, (0, 0, 1))]
        direction = [(30, 30, 30), (0, 0, 90)]
        self.assertTrue(len(crv) == len(direction))

        iges_model = Model()
        for k in range(len(crv)):
            sf = ExtrudedSurf(crv[k], direction[k])
            iges_model.clear()
            iges_model.add(sf.to_iges())
            iges_model.save('test_extrude-{}.igs'.format(k))

    def test_revolved(self):
        generatrix = [Line((50, 0, 0), (50, 0, 20)),
                      Arc.from_2pnt((100, 0, 0), (100, 0, 20), 180, (1, 0, 0))]
        center = [(20, 0, 0), (50, 0, 0)]
        axis = [(0, 0, 1), (0, 0, 1)]
        theta = [60, 180]

        iges_model = Model()
        for k in range(len(generatrix)):
            iges_model.clear()
            iges_model.add(generatrix[k].to_iges())
            sf = RevolvedSurf(center[k], axis[k], theta[k], generatrix[k])
            iges_model.add(sf.to_iges())
            iges_model.save('test_revolved-{}.igs'.format(k))
        self.assertTrue(True)

    def test_ruled(self):
        foil1 = WingProfile('M6', [(0, 0, 0), (10, 0, 0)]).crv
        foil2 = WingProfile('NACA0012', [(0, 0, 80), (15, 0, 80)]).crv
        surf = RuledSurf(foil1, foil2)
        model_file = Model()
        model_file.add(surf.to_iges())
        model_file.save('test_ruled.igs')
        self.assertTrue(True)

    def test_global_interp(self):
        # foil, N: num of section, L: length per section, p, q
        data = [('M6', 10, 0.7, 3, 3),
                ('M6', 10, 0.7, 3, 5),
                ('M6', 10, 0.7, 5, 3),
                ('M6', 10, 0.7, 5, 5),
                ('NACA0012', 10, 0.7, 3, 3),
                ('NACA0012', 10, 0.7, 3, 5),
                ('NACA0012', 10, 0.7, 5, 3),
                ('NACA0012', 10, 0.7, 5, 5),
                ('RAE2822', 10, 0.7, 3, 3),
                ('RAE2822', 10, 0.7, 3, 5),
                ('RAE2822', 10, 0.7, 5, 3),
                ('RAE2822', 10, 0.7, 5, 5)]

        iges_model = Model()
        for k in range(len(data)):
            foil = data[k][0]
            n = data[k][1]
            l = data[k][2]
            p = data[k][3]
            q = data[k][4]
            pts = read_airfoil_pts(foil)
            m, dim = pts.shape
            all_pts = np.zeros((n + 1, m, dim))
            for i in range(n + 1):
                all_pts[i] = np.copy(pts)
                for j in range(m):
                    all_pts[i][j][-1] = l * i
            wsf = GlobalInterpolatedSurf(all_pts, p, q)
            iges_model.clear()
            iges_model.add(wsf.to_iges())
            iges_model.save("test_global_interp-{}_{}_{}_{}.igs".format(k, foil, p, q))
        self.assertTrue(True)

    def test_local_interp(self):
        pass

    def test_coons(self):
        pass

    def test_skinned(self):
        # foil, N: num of section, L: length per section, p, q
        data = [('M6', 10, 0.7, 3, 3),
                ('M6', 10, 0.7, 3, 5),
                ('M6', 10, 0.7, 5, 3),
                ('M6', 10, 0.7, 5, 5),
                ('NACA0012', 10, 0.7, 3, 3),
                ('NACA0012', 10, 0.7, 3, 5),
                ('NACA0012', 10, 0.7, 5, 3),
                ('NACA0012', 10, 0.7, 5, 5),
                ('RAE2822', 10, 0.7, 3, 3),
                ('RAE2822', 10, 0.7, 3, 5),
                ('RAE2822', 10, 0.7, 5, 3),
                ('RAE2822', 10, 0.7, 5, 5)]

        iges_model = Model()
        for k in range(len(data)):
            foil = data[k][0]
            n = data[k][1]
            l = data[k][2]
            p = data[k][3]
            q = data[k][4]
            pts = read_airfoil_pts(foil)
            m, dim = pts.shape
            all_pts = np.zeros((n + 1, m, dim))
            for i in range(n + 1):
                all_pts[i] = np.copy(pts)
                for j in range(m):
                    all_pts[i][j][-1] = l * i

            crv_list = [GlobalInterpolatedCrv(all_pts[i], p) for i in range(n + 1)]
            wsf = Skinned(crv_list, p, q)
            iges_model.clear()
            iges_model.add(wsf.to_iges())
            iges_model.save("test_Skinned-{}_{}_{}_{}.igs".format(k, foil, p, q))
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
