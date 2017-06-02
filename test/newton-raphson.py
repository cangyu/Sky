import scipy
from scipy.optimize import newton
import math

A1 = 0.4
A2 = -1.2
A3 = 0.5


def p():
    return (1 - A1) * A3 * (scipy.exp(A2) - 1) / ((1 - A3) * A1 * A2 * scipy.exp(A2))


def f(x):
    return scipy.exp(x) - 1.0 - p() * x


def pf(x):
    return scipy.exp(x) - p()


pp = p()
p1z = math.log(pp)
fp1z = f(p1z)
A4 = newton(f, 10 * p1z, pf, maxiter=500)

print(pp, p1z, fp1z, A4, f(A4))
