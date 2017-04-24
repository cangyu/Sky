import numpy as np
import math
from src.msh.plot3d import *
from src.aircraft.wing import *

R = 50
U = 61
V = 15
W = 1


def calc_joint(xu_1: float, yu_1: float, xu: float, yu: float):
    r = np.array([xu - xu_1, yu - yu_1])
    a = r[0]
    b = r[1]
    t = a * xu + b * yu

    A = a * a + b * b
    B = -2 * t * a
    C = t * t - b * b * R * R
    delta = B * B - 4 * A * C

    x1 = (-B + math.sqrt(delta)) / (2 * A)
    x2 = (-B - math.sqrt(delta)) / (2 * A)
    y1 = (t - a * x1) / b
    y2 = (t - a * x2) / b

    return (x1, y1), (x2, y2)


def calc_msh_divide_pts(airfoil: Airfoil):
    ans = np.zeros((2, 2), float)

    cu = calc_joint(airfoil.x[-2], airfoil.y_up[-2], airfoil.x[-1], airfoil.y_up[-1])
    if (cu[0][1] > cu[1][1]):
        ans[0][0] = cu[0][0]
        ans[0][1] = cu[0][1]
    else:
        ans[0][0] = cu[1][0]
        ans[0][1] = cu[1][1]

    cd = calc_joint(airfoil.x[-2], airfoil.y_down[-2], airfoil.x[-1], airfoil.y_down[-1])
    if (cd[0][1] < cd[1][1]):
        ans[1][0] = cd[0][0]
        ans[1][1] = cd[0][1]
    else:
        ans[1][0] = cd[1][0]
        ans[1][1] = cd[1][1]

    return ans


def calc_msh_pts(airfoil: Airfoil):
    U = airfoil.pts_num
    pts = np.zeros((W, V, U, 3), float)

    # 起始坐标与角度
    itc = calc_msh_divide_pts(airfoil)
    uitc, ditc = itc[0], itc[1]

    sa = math.degrees(math.atan2(uitc[1], uitc[0]))
    if sa < 0:
        sa += 360

    ea = math.degrees(math.atan2(ditc[1], ditc[0]))
    if ea < 0:
        ea += 360

    ang_list = np.linspace(sa, ea, int(U))

    # 内层点
    for i in range(0, U):
        for j in range(0, 2):
            pts[0][0][i][j] = airfoil.pts[i][j]

    # 外层点
    pts[0][V - 1][0][0] = uitc[0]
    pts[0][V - 1][0][1] = uitc[1]
    pts[0][V - 1][U - 1][0] = ditc[0]
    pts[0][V - 1][U - 1][1] = ditc[1]
    for i in range(1, U - 1):
        ca = ang_list[i]
        pts[0][V - 1][i][0] = R * math.cos(math.radians(ca))
        pts[0][V - 1][i][1] = R * math.sin(math.radians(ca))

    # 内部点
    for i in range(0, U):
        for j in range(1, V - 1):
            ratio = j / (V - 1)
            for k in range(0, 2):
                pts[0][j][i][k] = pts[0][0][i][k] + (pts[0][V - 1][i][k] - pts[0][0][i][k]) * ratio

    # Output
    return Plot3D(U, V, W, pts)


if __name__ == '__main__':
    foil = 'NACA0012'
    naca0012 = Airfoil(foil)
    # calc_msh_pts(naca0012).output(foil)
    naca0012.iges()
