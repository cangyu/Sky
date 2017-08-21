import numpy as np
from random import random, randint
from collections.abc import Callable
from abc import abstractmethod


class Chromosome(object):
    def __init__(self, idx):
        """
        Individual representation used in Genetic Algorithm.
        :param idx: Index of this chromosome.
        :type idx: int
        """

        self.idx = idx
        self.fitness = float(0)
        self.value = float(0)

    @abstractmethod
    def calc_value(self, obj_func):
        pass

    @abstractmethod
    def calc_fitness(self, eval_func):
        pass


class RealCodedChromosome(Chromosome):
    def __init__(self, idx, n):
        """
        Parameters inside a chromosome are normalized to [0, 1].
        :param idx: Index of this chromosome.
        :type idx: int
        :param n: Number of parameters in this chromosome.
        :type n: int
        """

        super(RealCodedChromosome, self).__init__(idx)
        self.param = np.empty(n, float)
        for i in range(n):
            self.param[i] = random()

    def calc_value(self, obj_func):
        self.value = obj_func(self.param)

    def calc_fitness(self, eval_func):
        self.fitness = eval_func(self.param)

    def mutate(self):
        n = len(self.param)
        for i in range(n):
            self.param[i] = random()


class GA(object):
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
        self.individual = None


class RealCodedGA(GA):
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

        '''Parameter ranges'''
        self.param_range = np.empty((n, 2), float)

        for k, rg in enumerate(arg_rg):
            if len(rg) != 2 or rg[0] > rg[1]:
                raise AssertionError("Invalid input.")
            else:
                self.param_range[k] = rg

    @property
    def param_num(self):
        return len(self.param_range)

    def generate_population(self, n):
        """
        Generate a set of chromosome.
        :param n: Number of chromosome.
        :type n: int
        :return: List of chromosomes.
        """

        return [RealCodedChromosome(i, self.param_num) for i in range(n)]

    def find_optimal(self, n, rd, pc, pm):
        """
        Try to find the global optimal with given settings.
        :param n: Number of chromosome in the individual.
        :type n: int
        :param rd: Number of iteration.
        :type rd: int
        :param pc: Possibility of cross.
        :type pc: float
        :param pm: Possibility of mutate.
        :type pm: float
        :return: The global optimal chromosome.
        """

        '''Pre-check'''
        if n < 0 or rd < 0:
            raise AssertionError("Invalid input.")
        if pc < 0 or pc > 1:
            raise ValueError("Invalid \'pc\' setting.")
        if pm < 0 or pm > 1:
            raise ValueError("Invalid \'pm\' setting.")

        '''Init'''
        self.individual = self.generate_population(n)
        for k, chromosome in enumerate(self.individual):
            self.individual[k].calc_value(self.f_obj)
            self.individual[k].calc_fitness(self.f_eval)
        self.individual.sort()

        '''Iterate'''
        elite_num = 3
        n -= elite_num
        for generation in range(rd):
            next_gen = []

            '''Elitist Migration'''
            for i in range(elite_num):
                next_gen.append(self.individual[i])

            '''Tournament Selection'''
            for i in range(n):
                c1, c2 = 0, 0
                while c1 == c2:
                    c1 = randint(0, n)
                    c2 = randint(0, n)

                f1 = self.individual[c1].fitness
                f2 = self.individual[c2].fitness
                next_gen.append(self.individual[c1] if f1 > f2 else self.individual[c2])

            '''Cross-Over'''
            cross_flag = [False] * n
            for i in range(n):
                c1, c2 = 0, 0
                while c1 == c2 or cross_flag[c1] or cross_flag[c2]:
                    c1 = randint(0, n)
                    c2 = randint(0, n)

                cross_flag[c1] = cross_flag[c2] = True
                if random() < pc:
                    pos = randint(0, self.param_num)
                    next_gen[c1].param[pos], next_gen[c2].param[pos] = next_gen[c2].param[pos], next_gen[c1].param[pos]

            '''Mutation'''
            for i in range(n):
                if random() < pm:
                    next_gen[i].mutate()

            '''Sort for next generation'''
            self.individual = sorted(next_gen)
