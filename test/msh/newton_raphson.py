import scipy
from scipy.optimize import newton
import math

A1 = 0.4
A2 = -1.2
A3 = 0.5

p = (1 - A1) * A3 * (scipy.exp(A2) - 1) / ((1 - A3) * A1 * A2 * scipy.exp(A2))


def f(x):
    return scipy.exp(x) - 1.0 - p * x


def pf(x):
    return scipy.exp(x) - p


def ppf(x):
    return scipy.exp(x)


p1z = math.log(p)
fp1z = f(p1z)
A4 = newton(f, 5 * p1z, fprime=pf, maxiter=20, fprime2=ppf)

print(p, p1z, fp1z, A4, f(A4))
