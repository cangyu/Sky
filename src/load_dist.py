import math
from scipy.integrate import romberg
from abc import abstractmethod, ABCMeta

g = 9.80665


def q_inf(rho_inf, v_inf):
    return 0.5 * rho_inf * v_inf ** 2


def tons2newton(w):
    return w * 1000 * g


class LiftDist(object, metaclass=ABCMeta):
    def __init__(self, w, spn, rho, v, root_lift):
        self.rho_inf = rho
        self.v_inf = v
        self.q_inf = q_inf(rho, v)
        self.m_inf = rho * v
        self.span = spn
        self.span2 = spn / 2
        self.payload = tons2newton(w)
        self.root_lift = root_lift

    @abstractmethod
    def lift_at(self, u):
        """
        Y=rho*v*gamma
        :param u: Relative position.
        :return: 2D lift at given position.
        :rtype: float
        """

        pass

    def gamma_at(self, rel_pos):
        return self.lift_at(rel_pos) / self.m_inf

    def cl_at(self, rel_pos, b):
        return self.lift_at(rel_pos) / (self.q_inf * b)

    def lift_between(self, rel_begin, rel_end):
        return romberg(self.lift_at, rel_begin, rel_end) * self.span2


class EllipticLiftDist(LiftDist):
    def __init__(self, w, spn, rho, v):
        """
        Elliptic Lift distribution in span-wise direction.
        This is the ideal distribution for reducing induced drag.
        :param w: Weight of the aircraft in cruise, in tons.
        :type w: float
        :param spn: Span of the aircraft, in meters.
        :type spn: float
        :param rho: Density of the free-stream.
        :type rho: float
        :param v: Velocity of the free-stream.
        :type v: float
        """

        root_lift = 4 * tons2newton(w) / (math.pi * spn)
        super(EllipticLiftDist, self).__init__(w, spn, rho, v, root_lift)

    def lift_at(self, u):
        return self.root_lift * math.sqrt(1 - u ** 2)


class LinearLiftDist(LiftDist):
    def __init__(self, w, spn, rho, v):
        """
        Triangular/Linear Lift distribution in span-wise direction.
        :param w: Weight of the aircraft in cruise, in tons.
        :type w: float
        :param spn: Span of the aircraft, in meters.
        :type spn: float
        :param rho: Density of the free-stream.
        :type rho: float
        :param v: Velocity of the free-stream.
        :type v: float
        """

        root_lift = 2 * tons2newton(w) / spn
        super(LinearLiftDist, self).__init__(w, spn, rho, v, root_lift)

    def lift_at(self, u):
        return self.root_lift * (1 - u)


class UniformLiftDist(LiftDist):
    def __init__(self, w, spn, rho, v):
        """
        Rectangular/Uniform Lift distribution in span-wise direction.
        :param w: Weight of the aircraft in cruise, in tons.
        :type w: float
        :param spn: Span of the aircraft, in meters.
        :type spn: float
        :param rho: Density of the free-stream.
        :type rho: float
        :param v: Velocity of the free-stream.
        :type v: float
        """

        root_lift = tons2newton(w) / spn
        super(UniformLiftDist, self).__init__(w, spn, rho, v, root_lift)

    def lift_at(self, u):
        return self.root_lift
