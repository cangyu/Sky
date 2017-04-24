import numpy as np
import math
from src.msh.plot3d import *

R_MIN = 50
R_MAX = 100
THETA_MIN = 30
THETA_MAX = 120
U = 61
V = 16
W = 20

fn = "../../result/cylinder.xyz"

if __name__ == "__main__":

    pts = np.zeros((W, V, U, 3), float)
    for w in range(0, W):
        cw = w
        for v in range(0, V):
            cr = R_MIN + v / (V - 1) * (R_MAX - R_MIN)
            for u in range(0, U):
                ct = math.radians(THETA_MIN + u / (U - 1) * (THETA_MAX - THETA_MIN))
                pts[w][v][u][0] = cr * math.cos(ct)
                pts[w][v][u][1] = cr * math.sin(ct)
                pts[w][v][u][2] = cw

    msh = Plot3D(U, V, W, pts)
    msh.output(fn)
