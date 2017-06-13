import numpy as np
from src.msh.tfi import Linear_TFI_2D
from src.msh.fluent import XF_MSH


def rectangular(L: float, W: float):
    return Linear_TFI_2D(lambda u: np.array([L * u, 0, 0]),
                         lambda v: np.array([0, W * v, 0]),
                         lambda u: np.array([L * u, W, 0]),
                         lambda v: np.array([L, W * v, 0]))


U = 2
V = 2
L = 5
W = 2

if __name__ == '__main__':
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    ppu, ppv = np.meshgrid(u_list, v_list, indexing='xy')
    msh = rectangular(L, W)
    msh.calc_grid(ppu, ppv)
    fmsh = XF_MSH.from_str2d(msh.get_grid())
    fmsh.save('test.msh')
