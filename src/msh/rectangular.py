import numpy as np
import math
from src.msh.plot3d import *

X_MIN = 0
X_MAX = 100
Y_MIN = 0
Y_MAX = 60
Z_MIN = 0
Z_MAX = 40

U = 61
V = 16
W = 20

fn = "../../result/rectangular.xyz"


def p2r(i, I, min, max):
    return min + (max - min) * i / (I - 1)


if __name__ == "__main__":

    pts = np.zeros((W, V, U, 3), float)
    for w in range(0, W):
        cw = p2r(w, W, Z_MIN, Z_MAX)
        for v in range(0, V):
            cv = p2r(v, V, Y_MIN, Y_MAX)
            for u in range(0, U):
                cu = p2r(u, U, X_MIN, X_MAX)
                pts[w][v][u][0] = cu
                pts[w][v][u][1] = cv
                pts[w][v][u][2] = cw

    msh = Plot3D(U, V, W, pts)
    msh.output(fn)
