import math
import numpy as np
import scipy
from scipy.optimize import newton


def single_exponential(begin, end, N, A):
    rho = np.linspace(0, 1.0, N)
    r = np.copy(rho)
    for i in range(len(r)):
        r[i] = (math.exp(A * r[i]) - 1) / (math.exp(A) - 1)

    delta = end - begin
    for i in range(len(r)):
        r[i] *= delta
        r[i] += begin

    return r


def double_exponential(begin, end, N, A1, A2, A3):
    rho = np.linspace(0, 1.0, N)
    p = (1 - A1) * A3 * (scipy.exp(A2) - 1) / ((1 - A3) * A1 * A2 * scipy.exp(A2))
    if p <= 0:
        raise ValueError("Invalid parameters!")

    def f(x):
        return scipy.exp(x) - 1.0 - p * x

    def pf(x):
        return scipy.exp(x) - p

    p1z = math.log(p)
    A4 = newton(f, 10 * p1z, pf, maxiter=500)

    def r(x):
        if x <= A3:
            return A1 * (math.exp(A2 / A3 * x) - 1) / (math.exp(A2) - 1)
        else:
            return A1 + (1 - A1) / (math.exp(A4) - 1) * (math.exp(A4 / (1 - A3) * (x - A3)) - 1)

    ans = np.copy(rho)
    for i in range(len(ans)):
        ans[i] = r(ans[i])

    delta = end - begin
    for i in range(len(ans)):
        ans[i] *= delta
        ans[i] += begin

    return ans
