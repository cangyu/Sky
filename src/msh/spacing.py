import math
import numpy as np


def single_exponential(begin, end, N, A):
    rho = np.linspace(begin, end, N + 1)
    r = np.copy(rho)
    for i in range(len(r)):
        r[i] = (math.exp(A * r[i]) - 1) / (math.exp(A) - 1)
    return r
