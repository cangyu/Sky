import numpy as np


def normalize(vector):
    tmp = np.sqrt(sum(map(lambda a: a ** 2, vector)))
    return vector / tmp


class DCM(object):
    def __init__(self, base1, base2):
        """
        方向余弦矩阵
        :param base1: 起始坐标轴标架 
        :param base2: 目标坐标轴标架
        """

        I = normalize(base1[0])
        J = normalize(base1[1])
        K = normalize(base1[2])
        i = normalize(base2[0])
        j = normalize(base2[1])
        k = normalize(base2[2])

        self.dcm = np.matrix([[np.dot(I, i), np.dot(I, j), np.dot(I, k)],
                              [np.dot(J, i), np.dot(J, j), np.dot(J, k)],
                              [np.dot(K, i), np.dot(K, j), np.dot(K, k)]])

    @property
    def rot_matrix(self):
        return self.dcm
