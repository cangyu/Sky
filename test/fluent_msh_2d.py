import numpy as np
from src.msh.tfi import Linear_TFI_2D
from src.msh.fluent import XF_MSH, BCType


def rectangular(U: int, V: int, L: float, W: float):
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    ppu, ppv = np.meshgrid(u_list, v_list, indexing='ij')
    msh = Linear_TFI_2D(lambda u: np.array([L * u, 0, 0]),
                        lambda v: np.array([0, W * v, 0]),
                        lambda u: np.array([L * u, W, 0]),
                        lambda v: np.array([L, W * v, 0]))
    msh.calc_grid(ppu, ppv)
    fluent_grid = XF_MSH.from_str2d(msh.get_grid())
    fluent_grid.save('rect_{}_{}.msh'.format(U, V))


def curve_rect(U: int, V: int, L: float, H1: float, H2: float, H3: float):
    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    ppu, ppv = np.meshgrid(u_list, v_list, indexing='ij')
    msh = Linear_TFI_2D(lambda u: np.array([u * L, 4 * H3 * u * (1 - u), 0]),
                        lambda v: np.array([0, v * H1, 0]),
                        lambda u: np.array([u * L, (H1 * (1 - u * u) + H2 * u * u), 0]),
                        lambda v: np.array([L, v * H2, 0]))
    msh.calc_grid(ppu, ppv)
    fluent_grid = XF_MSH.from_str2d(msh.get_grid())
    fluent_grid.save('crv_rect_{}_{}.msh'.format(U, V))


def airfoil():
    pass


if __name__ == '__main__':
    rectangular(30, 10, 3, 1)
    curve_rect(50, 25, 100, 40, 60, 10)
