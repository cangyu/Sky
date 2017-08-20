import numpy as np


class LatinHyperCube(object):
    def __init__(self, entry_list):
        """
        Latin hyper cube sampling.
        :param entry_list: Parameter entries.
        """

        '''Parameter entry list'''
        self.entry_list = entry_list

        '''Defensive check'''
        if len(self.entry_list) == 0:
            raise ValueError("Invalid input.")

        '''Random permutation'''
        self.rp = None
        self.n = len(self.entry_list)
        self.m = len(self.entry_list[0])

        '''Defensive check'''
        for k, sp in enumerate(self.entry_list):
            if len(sp) != self.m:
                raise AssertionError("Inconsistent length of entry {}.".format(k))

        '''Display distribution'''
        self.r1, self.r2 = np.meshgrid(np.arange(self.m), np.arange(self.n))
        self.data = np.zeros((self.n, self.m), int)

    def sample(self):
        """
        Sampling on given parameters.
        :return: Random sampling.
        """

        self.rp = np.empty((self.n, self.m), int)
        for i in range(self.n):
            self.rp[i] = np.random.permutation(self.m)

        ret = []
        for i in range(self.m):
            cur_comb = []
            for j in range(self.n):
                cur_comb.append(self.entry_list[j][self.rp[j][i]])
            ret.append(cur_comb)

        return ret
