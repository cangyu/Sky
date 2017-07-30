import math
from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate


class WingFrame(object):
    def __init__(self, xf, xt, yf, yt, z):
        """
        机翼外框描述
        :param xf: 前缘x坐标的参数方程
        :param xt: 后缘x坐标的参数方程
        :param yf: 前缘y坐标的参数方程
        :param yt: 后缘y坐标的参数方程
        :param z: z坐标的参数方程
        """

        self.f_xf = deepcopy(xf)
        self.f_xt = deepcopy(xt)
        self.f_yf = deepcopy(yf)
        self.f_yt = deepcopy(yt)
        self.f_z = deepcopy(z)

    def x_front(self, u):
        return self.f_xf(u)

    def x_tail(self, u):
        return self.f_xt(u)

    def y_front(self, u):
        return self.f_yf(u)

    def y_tail(self, u):
        return self.f_yt(u)

    def z(self, u):
        return self.f_z(u)

    def show(self, n=1000):
        u_dist = np.linspace(0, 1.0, n)
        z = np.empty(n, float)
        xf = np.empty(n, float)
        xt = np.empty(n, float)
        for k in range(n):
            z[k] = self.z(u_dist[k])
            xf[k] = self.x_front(u_dist[k])
            xt[k] = self.x_tail(u_dist[k])

        plt.figure()
        plt.plot(z, xf, label='Front')
        plt.plot(z, xt, label='Tail')
        plt.legend()
        plt.gca().invert_yaxis()
        plt.show()


class BWBFrame(WingFrame):
    def __init__(self, c_root, c_mid, c_tip, b_mid, b_tip, alpha_root, alpha_mid, alpha_tip):
        """
        BWB(Blended Wing Body)构型飞行器参数化描述
        :param c_root: 根部弦长
        :type c_root: float
        :param c_mid: 中间弦长
        :type c_mid: float
        :param c_tip: 翼尖弦长
        :type c_tip: float
        :param b_mid: 内翼宽度
        :type b_mid: float
        :param b_tip: 机翼半展长
        :type b_tip: float
        :param alpha_root: 内翼平均后掠角
        :type alpha_root: float
        :param alpha_mid: 中段平均后掠角
        :type alpha_mid: float
        :param alpha_tip: 翼尖平均后掠角
        :type alpha_tip: float
        """

        self.Cr = c_root
        self.Cm = c_mid
        self.Ct = c_tip
        self.Bm = b_mid
        self.Bt = b_tip
        self.Dm = d_mid
        self.Dt = d_tip
        self.Ar = alpha_root
        self.Am = alpha_mid
        self.At = alpha_tip

        '''Calculate pivots on each curve'''
        front_pnt = np.empty((3, 3), float)
        tail_pnt = np.empty((3, 3), float)
        front_pnt[0] = np.zeros(3)
        tail_pnt[0] = np.array([self.Cr, 0, 0])
        front_pnt[1] = np.array([self.Bm * math.tan(math.radians(self.Am)), 0, self.Bm])
        tail_pnt[1] = np.array([front_pnt[1][0] + self.Cm, 0, self.Bm])
        front_pnt[2] = np.array([front_pnt[1][0] + (self.Bt - self.Bm) * math.tan(math.radians(self.At)), 0, self.Bt])
        tail_pnt[2] = np.array([front_pnt[2][0] + self.Ct, 0, self.Bt])

        '''Build interpolated functions'''
        u = np.array([0, self.Bm / self.Bt, 1.0])
        z = np.array([0, self.Bm, self.Bt])
        xf = np.array([front_pnt[0][0], front_pnt[1][0], front_pnt[2][0]])
        yf = np.array([front_pnt[0][1], front_pnt[1][1], front_pnt[2][1]])
        xt = np.array([tail_pnt[0][0], tail_pnt[1][0], tail_pnt[2][0]])
        yt = np.array([tail_pnt[0][1], tail_pnt[1][1], tail_pnt[2][1]])

        super(BWBFrame, self).__init__(interpolate.make_interp_spline(u, xf, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       interpolate.make_interp_spline(u, xt, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       interpolate.make_interp_spline(u, yf, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       interpolate.make_interp_spline(u, yt, 3, bc_type=([(1, 0)], [(2, 0)])),
                                       lambda t: z[1] * t / u[1] if t <= u[1] else z[1] + (z[2] - z[1]) * (t - u[1]) / (u[2] - u[1]))


def chebshev_dist(start, end, n):
    """
    生成切比雪夫点
    :param start: 起始值
    :type start: float
    :param end: 终止值
    :type end: float
    :param n: 采样点数量
    :type n: int
    :return: n个点
    """

    ang = np.linspace(math.pi, 0, n)
    pr = np.zeros(n)
    for i in range(0, n):
        pr[i] = math.cos(ang[i])
        pr[i] = start + (end - start) / 2 * (pr[i] + 1)

    return pr
