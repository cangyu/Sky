import unittest
import sys
import math
import numpy as np
from numpy.linalg import inv, det
from random import random, randint
from collections.abc import Callable
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
Latin Hypercube Sampling
"""


class LHS(object):
    def __init__(self, entry_list):
        """
        Latin Hypercube Sampling.
        :param entry_list: Parameter entries.
        """

        '''Defensive check'''
        if len(entry_list) == 0:
            raise ValueError("Invalid input.")
        n = len(entry_list[0])
        for k in range(1, len(entry_list)):
            assert len(entry_list[k]) == n

        '''Parameter entry list'''
        self.entry = np.copy(entry_list)

    @property
    def var_num(self):
        return len(self.entry)

    @property
    def sample_num(self):
        return len(self.entry[0])

    def sample(self):
        """
        Sampling on given parameters.
        :return: Random sampling.
        """

        rp = np.copy([np.random.permutation(self.sample_num) for i in range(self.var_num)]).transpose()
        return [[self.entry[j][rp[i][j]] for j in range(self.var_num)] for i in range(self.sample_num)]


class Person(object):
    def __init__(self, _nm, _ag):
        self._name = _nm
        self._age = _ag

    @property
    def name(self):
        return self._name

    @property
    def age(self):
        return self._age

    def __repr__(self):
        return "My name is {}, I'm {} years old now!".format(self.name, self.age)


class LHSTester(unittest.TestCase):
    def test_sample(self):
        a = (0, 1, 2, 3, 4, 5)
        b = ['a', 'b', 'c', 'd', 'e', 'f']
        c = ('I', 'II', 'III', 'IV', 'V', 'VI')
        d = np.array([1.1, 2.2, 3.3, 4.4, 5.5, 6.6])
        e = (Person('ggsmd', 21), Person('sbtty', 22), Person('smdhx', 23), Person('shzyh', 24), Person('sgdb', 25), Person('tltbpa', 26))
        self.assertTrue(len(a) == len(b) == len(c) == len(d) == len(e))

        lhc = LHS([a, b, c, d, e])
        sp = lhc.sample()
        for i in range(len(sp)):
            print('\n {}'.format(sp[i]))


"""
Kriging Surrogate Model
"""


class Kriging(object):
    def __init__(self, x, z):
        """
        Kriging Surrogate Model.
        :param x: Sample points.
        :param z: Values on sample points.
        """

        '''Pre-Check'''
        n = len(x)
        if len(x) != len(z) or n == 0:
            raise AssertionError("Invalid input.")
        nd = len(x[0])
        f = np.ones(n, float)

        super(Kriging, self).__init__()

        '''Sample points and values'''
        self.x = np.copy(x)
        self.z = np.copy(z)

        '''Distance between sample points'''
        d = np.empty((n, n, nd), float)
        for i in range(n):
            for j in range(i, n):
                for k in range(nd):
                    d[i][j][k] = d[j][i][k] = math.fabs(self.x[i][k] - self.x[j][k])

        def mle(_t):
            _r = np.empty((n, n), float)
            for _ci in range(n):
                for _cj in range(_ci, n):
                    _pw = 0
                    for _ck in range(nd):
                        _pw += _t[_ck] * d[_ci][_cj][_ck] ** 2
                    _r[_ci][_cj] = _r[_cj][_ci] = math.exp(-_pw)

            _r_inv = inv(_r)
            _tmp1 = np.dot(f, _r_inv)
            _expect = np.dot(_tmp1, self.z) / np.dot(_tmp1, f)
            _tmp2 = self.z - _expect * f
            _variance = np.dot(np.dot(_tmp2, _r_inv), _tmp2) / n
            return -(n * math.log(_variance) + math.log(det(_r))) / 2

        def pmle(_t):
            _param = np.copy(_t)
            for _k, _p in enumerate(_param):
                _param[_k] = math.pow(10, _p)
            return mle(_param)

        '''Determine theta_k using GA'''
        rg = np.copy([(-3, 1.3)] * nd)
        rga = RealCodedGA(rg, pmle, pmle)
        self.theta = rga.find_optimal(100 * nd, 20 * nd)
        for _k, _p in enumerate(self.theta):
            self.theta[_k] = math.pow(10, _p)

        '''Correlation under Gauss function'''
        self.R = np.empty((n, n), float)
        for i in range(n):
            for j in range(i, n):
                pdx = 0
                for k in range(nd):
                    pdx += self.theta[k] * d[i][j][k] ** 2
                self.R[i][j] = self.R[j][i] = math.exp(-pdx)

        '''Expect and variance'''
        r_inv = inv(self.R)
        tmp1 = np.dot(f, r_inv)
        e = np.dot(tmp1, self.z) / np.dot(tmp1, f)
        tmp2 = self.z - e * f
        self.var = np.inner(np.dot(tmp2, r_inv), tmp2) / n

        '''Covariance'''
        self.cov = self.var * self.R

        '''Semi-Variance'''
        self.r = np.full((n, n), self.var) - self.cov

    def interp(self, x0):
        """
        Calculate the response value at given point.
        :param x0: Observation point.
        :return: Response value.
        """

        n = len(self.x)
        nd = len(x0)

        r0 = np.empty(n, float)
        for i in range(n):
            tmp = 0.
            for d in range(nd):
                tmp += self.theta[d] * (x0[d] - self.x[i][d]) ** 2
            r0[i] = self.var - self.var * math.exp(-tmp)

        '''Matrix Coefficients'''
        A = np.empty((n + 1, n + 1), float)
        for i in range(n):
            for j in range(n):
                A[i][j] = self.r[i][j]
        for i in range(n):
            A[-1][i] = 1
        for j in range(n):
            A[j][-1] = 1
        A[-1][-1] = 0

        '''RHS'''
        b = np.empty(n + 1, float)
        for i in range(n):
            b[i] = r0[i]
        b[-1] = 1

        '''Solve linear system: 'Ax = b' '''
        x = np.dot(b, inv(A.transpose()))

        '''Assemble'''
        ans = 0
        for i in range(n):
            ans += x[i] * self.z[i]

        return ans

    def show(self):
        pass


class Branin(object):
    def __init__(self, a=1.0, b=5.1 / (4 * np.pi ** 2), c=5 / np.pi, r=6.0, s=10.0, t=1 / (8 * np.pi)):
        """
        Branin-Hoo function which has 3 global minima.
        """
        self.a = a
        self.b = b
        self.c = c
        self.r = r
        self.s = s
        self.t = t

    def __call__(self, x):
        """
        Calculate the value at (x1, x2).
        :param x: Parameter vector.
        :return: Value at (x1, x2).
        :rtype: float
        """

        x1 = x[0]
        x2 = x[1]
        ans = self.s
        ans += self.s * (1.0 - self.t) * math.cos(x1)
        ans += self.a * (x2 - self.b * x1 ** 2 + self.c * x1 - self.r) ** 2
        return ans


class KrigingTester(unittest.TestCase):
    def test_interp(self):
        branin_func = Branin()
        lhc = LHS(np.array([np.linspace(-5, 10, 20), np.linspace(0, 15, 20)]))

        x = lhc.sample()
        y = np.empty(len(x), float)
        for k, vx in enumerate(x):
            y[k] = branin_func(vx)

        kg = Kriging(x, y)
        x0 = (-math.pi, 12.275)
        x1 = (math.pi, 2.275)
        x2 = (9.42478, 2.475)
        print("\nGlobal minimum: {} Kriging: {} Actual: {}".format(x0, kg.interp(x0), branin_func(x0)))
        print("\nGlobal minimum: {} Kriging: {} Actual: {}".format(x1, kg.interp(x1), branin_func(x1)))
        print("\nGlobal minimum: {} Kriging: {} Actual: {}".format(x2, kg.interp(x2), branin_func(x2)))
        self.assertTrue(np.allclose(list(map(lambda u: kg.interp(u), x)), y))


"""
Genetic Algorithm
"""


class Chromosome(object):
    def __init__(self):
        """
        Individual representation used in Genetic Algorithm.
        """

        self.param = None
        self.fitness = float(0)
        self.value = float(0)

    def __lt__(self, other):
        return self.fitness > other.fitness


class RealCodedChromosome(Chromosome):
    def __init__(self, n):
        """
        Parameters inside a chromosome are normalized to [0, 1].
        :param n: Number of parameters in this chromosome.
        :type n: int
        """

        super(RealCodedChromosome, self).__init__()
        self.param = np.empty(n, float)
        for i in range(n):
            self.param[i] = random()


class GeneticAlgorithm(object):
    def __init__(self, obj_func, eval_func):
        """
        Genetic Algorithm.
        :param obj_func: The object function.
        :type obj_func: Callable
        :param eval_func: The evaluation function.
        :type eval_func: Callable
        """

        '''Object function'''
        self.f_obj = obj_func

        '''Evaluation function'''
        self.f_eval = eval_func

        '''Population set'''
        self.cur_generation = []


class RealCodedGA(GeneticAlgorithm):
    def __init__(self, arg_rg, obj_func, eval_func):
        """
        Real coded Genetic Algorithm.
        :param arg_rg: List of parameter ranges.
        :param obj_func: The object function.
        :type obj_func: Callable
        :param eval_func: The evaluation function.
        :type eval_func: Callable
        """

        super(RealCodedGA, self).__init__(obj_func, eval_func)

        '''Defensive check'''
        n = len(arg_rg)
        if n == 0:
            raise ValueError("Invalid input.")
        self.param_num = n

        '''Parameter ranges'''
        self.param_range = np.empty((n, 2), float)
        for k, rg in enumerate(arg_rg):
            if len(rg) != 2:
                raise AssertionError("Invalid input.")
            self.param_range[k] = rg

    def param_transform(self, p):
        ret = np.empty(self.param_num)
        for k in range(self.param_num):
            ret[k] = self.param_range[k][0] + (self.param_range[k][1] - self.param_range[k][0]) * p[k]
        return ret

    def find_optimal(self, n, rd, pm=0.05):
        """
        Try to find the global optimal with given settings.
        :param n: Number of chromosome in the cur_generation.
        :type n: int
        :param rd: Number of iteration.
        :type rd: int
        :param pm: Possibility of mutate.
        :type pm: float
        :return: The global optimal chromosome.
        """

        '''Regard top 1% individuals as elites'''
        elite_num = int(n / 100)

        '''Pre-check'''
        if n < elite_num or rd < 0:
            raise AssertionError("Invalid input.")
        if pm < 0 or pm > 1:
            raise ValueError("Invalid \'pm\' setting.")

        gen_digits = int(math.log10(rd)) + 1

        def report(_gen):
            if sys.flags.debug:
                _cur_best = self.cur_generation[0]
                _param = self.param_transform(_cur_best.param)
                _val = _cur_best.value
                print("Generation {0:{1}}: Best individual: {2} Value: {3}".format(_gen, gen_digits, _param, _val))

        '''Init'''
        self.cur_generation = []
        for i in range(n):
            c = RealCodedChromosome(self.param_num)
            real_param = self.param_transform(c.param)
            c.value = self.f_obj(real_param)
            c.fitness = self.f_eval(real_param)
            self.cur_generation.append(c)
        self.cur_generation.sort()
        report(0)

        def select(pop, lo, hi):
            """
            Tournament Selection.
            :return: Random selected chromosome.
            """

            c1, c2 = 0, 0
            while c1 == c2:
                c1 = randint(lo, hi)
                c2 = randint(lo, hi)

            return pop[c1] if pop[c1].fitness > pop[c2].fitness else pop[c2]

        def cross(p1, p2):
            ratio = random()
            ret = RealCodedChromosome(self.param_num)
            ret.param = ratio * p1.param + (1 - ratio) * p2.param
            return ret

        def mutate(idv):
            for j in range(self.param_num):
                if random() < pm:
                    idv.param[j] = random()

        '''Iterate'''
        for generation in range(1, rd + 1):
            next_gen = []

            '''Elite migration'''
            for i in range(elite_num):
                next_gen.append(self.cur_generation[i])

            '''Genetic operations'''
            for i in range(n - elite_num):
                id1 = select(self.cur_generation, 0, n - 1)
                id2 = select(self.cur_generation, 0, n - 1)
                offspring = cross(id1, id2)
                mutate(offspring)
                next_gen.append(offspring)

            '''Update and sort'''
            for k, chromosome in enumerate(next_gen):
                real_param = self.param_transform(chromosome.param)
                next_gen[k].value = self.f_obj(real_param)
                next_gen[k].fitness = self.f_eval(real_param)
            self.cur_generation = sorted(next_gen)
            report(generation)

        return self.param_transform(self.cur_generation[0].param)


class NashRealCodedGA(RealCodedGA):
    def __init__(self, arg_rg, obj_func, eval_func):
        """
        Nash Real coded Genetic Algorithm.
        :param arg_rg: List of parameter ranges.
        :param obj_func: The object function.
        :type obj_func: Callable
        :param eval_func: The evaluation function.
        :type eval_func: Callable
        """

        super(NashRealCodedGA, self).__init__(arg_rg, obj_func, eval_func)


class Rastrigin(object):
    def __call__(self, x):
        n = len(x)
        a = 10
        return a * n + sum(list(map(lambda u: u ** 2 - a * math.cos(2 * math.pi * u), x)))

    def show(self):
        x1 = np.linspace(-5.12, 5.12, 100)
        x2 = np.linspace(-5.12, 5.12, 100)

        X, Y = np.meshgrid(x1, x2, indexing='ij')

        val = np.empty((100, 100))
        for i in range(100):
            for j in range(100):
                val[i][j] = self.__call__((x1[i], x2[j]))

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot_surface(X, Y, val, rstride=1, cstride=1, cmap=plt.cm.hot)
        plt.show()


class Ackley(object):
    def __call__(self, x):
        return -20 * math.exp(-0.2 * math.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))) - math.exp(0.5 * (math.cos(2 * math.pi * x[0]) + math.cos(2 * math.pi * x[1]))) + math.e + 20

    def show(self):
        x1 = np.linspace(-5.12, 5.12, 100)
        x2 = np.linspace(-5.12, 5.12, 100)

        gx, gy = np.meshgrid(x1, x2, indexing='ij')

        val = np.empty((100, 100))
        for i in range(100):
            for j in range(100):
                val[i][j] = self.__call__((x1[i], x2[j]))

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot_surface(gx, gy, val, rstride=1, cstride=1, cmap=plt.cm.hot)
        plt.show()


class Sphere(object):
    def __call__(self, x):
        return sum(list(map(lambda u: u ** 2, x)))


class Rosenbrock(object):
    def __call__(self, x):
        n = len(x)
        if n < 2:
            raise AssertionError("Insufficient input parameters.")

        ret = 0
        for i in range(n - 1):
            ret += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (x[i] - 1) ** 2

        return ret


class GATester(unittest.TestCase):
    def test_rastrigin(self):
        rastrigin = Rastrigin()
        f_rastrigin = rastrigin.__call__
        rg = np.array([[-5.12, 5.12],
                       [-5.12, 5.12],
                       [-5.12, 5.12],
                       [-5.12, 5.12]])
        print('Testing Rastrigin function with {} variables ...'.format(len(rg)))
        rga = RealCodedGA(rg, f_rastrigin, lambda u: -f_rastrigin(u))
        ans = rga.find_optimal(600, 100, 0.05)
        print('Global Minimum: {}, Param: {}\n'.format(f_rastrigin(ans), ans))
        self.assertTrue(True)

    def test_ackley(self):
        ackley = Ackley()
        f_ackley = ackley.__call__
        rg = np.array([[-5.12, 5.12],
                       [-5.12, 5.12]])
        print('Testing Ackley function ...')
        rga = RealCodedGA(rg, f_ackley, lambda u: -f_ackley(u))
        ans = rga.find_optimal(300, 60, 0.05)
        print('Global Minimum: {}, Param: {}\n'.format(f_ackley(ans), ans))
        self.assertTrue(True)

    def test_sphere(self):
        sphere = Sphere()
        f_sphere = sphere.__call__
        rg = np.array([[-5.12, 5.12],
                       [-5.12, 5.12]])
        print('Testing Sphere function with {} variables ...'.format(len(rg)))
        rga = RealCodedGA(rg, f_sphere, lambda u: -f_sphere(u))
        ans = rga.find_optimal(300, 60, 0.05)
        print('Global Minimum: {}, Param: {}\n'.format(f_sphere(ans), ans))
        self.assertTrue(True)

    def test_rosenbrock(self):
        rosenbrock = Rosenbrock()
        f_rosenbrock = rosenbrock.__call__
        rg = np.array([[-5.12, 5.12],
                       [-5.12, 5.12],
                       [-5.12, 5.12]])
        print('Testing Rosenbrock function with {} variables ...'.format(len(rg)))
        rga = RealCodedGA(rg, f_rosenbrock, lambda u: -f_rosenbrock(u))
        ans = rga.find_optimal(1000, 300, 0.42)
        print('Global Minimum: {}, Param: {}\n'.format(f_rosenbrock(ans), ans))
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
