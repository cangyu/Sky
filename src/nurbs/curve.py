import numpy as np
import scipy
import math
from src.nurbs.point import Point2D

class Curve2D(object):
    def __init__(self, pts, p=5):
        """
        p阶，n+1个控制点，m+1个节点,全局插值
        """
        self.pts = pts
        self.p = p
        self.n = len(self.pts) - 1
        self.m = self.n + self.p + 1
        self.param = np.zeros(self.n + 1, float)
        self.knots = np.zeros(self.m + 1, float)
        self.coef = np.zeros((self.n + 1, self.n + 1), float)

    def calc_param(self, method='centripetal'):
        self.param[0] = 0.0
        self.param[self.n] = 1.0

        if method not in ['uniform', 'chord', 'centripetal']:
            raise Exception("Invalid parameter!")

        if method == 'uniform':
            inc = 1 / self.n
            for i in range(1, self.n):
                self.param[i] = self.param[i - 1] + inc
        else:
            dist = np.zeros(self.n + 1, float)
            for i in range(1, self.n + 1):
                dist[i] = Point2D.dist(self.pts[i], self.pts[i - 1])

            d = 0
            if method == 'chord':
                for i in range(1, self.n + 1):
                    d += dist[i]
            else:
                for i in range(1, self.n + 1):
                    dist[i] = math.sqrt(dist[i])
                    d += dist[i]

            for i in range(1, self.n):
                self.param[i] = self.param[i - 1] + dist[i] / d

    def calc_knots(self):
        """
        取平均值方法计算节点
        """
        for i in range(0, self.p + 1):
            self.knots[i] = 0.0
            self.knots[self.m - i] = 1.0

        acc = 0.0
        for i in range(0, self.p):
            acc += self.param[i]

        for j in range(1, self.n - self.p):
            acc -= self.param[j - 1]
            acc += self.param[self.p - 1 + j]
            self.knots[self.p + j] = acc / self.p

    def calc_coef(self):
        """
        计算系数矩阵
        """
        pass

