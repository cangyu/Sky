import numpy as np
from math import log10
from random import random, randint
from collections.abc import Callable


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

    def find_optimal(self, n, rd, pm):
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

        gen_digits = int(log10(rd)) + 1

        def report(_gen):
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

        return self.cur_generation[0]


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

    def find_optimal(self, n, rd, pm):
        pass
