import math
import numpy as np
import matplotlib.pyplot as plt

N = 100
rho = np.linspace(0, 1.0, N + 1)

A = 2
r = np.copy(rho)
for i in range(len(r)):
    r[i] = (math.exp(A * r[i]) - 1) / (math.exp(A) - 1)

plt.figure()
plt.plot(rho, r)
plt.show()
