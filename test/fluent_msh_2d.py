import numpy as np
from src.msh.tfi import Linear_TFI_2D
from src.msh.fluent import XF_MSH, BCType
from src.nurbs.surface import Coons
from src.nurbs.curve import Arc, Line
from src.aircraft.wing import WingProfile


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


def airfoil(U, V, foil, ends, thk, order, arc_start, arc_end, theta, nv):
    u0 = WingProfile(foil, ends, thk, order).nurbs_rep
    u1 = Arc.from_2pnt(arc_start, arc_end, theta, nv)
    v0 = Line(u0.start, u1.start)
    v1 = Line(u0.end, u1.end)

    u_list = np.linspace(0, 1.0, U + 1)
    v_list = np.linspace(0, 1.0, V + 1)
    ppu, ppv = np.meshgrid(u_list, v_list, indexing='ij')
    msh = Linear_TFI_2D(u0, v0, u1, v1)
    msh.calc_grid(ppu, ppv)
    fluent_grid = XF_MSH.from_str2d(msh.get_grid())
    fn = '{}_{}_{}.msh'.format(foil, U, V)
    fluent_grid.save(fn)


if __name__ == '__main__':
    rectangular(30, 10, 3, 1)
    curve_rect(50, 25, 100, 40, 60, 10)
    airfoil(130, 45, 'M6', [[0, 0, 5], [20, 0, 5]], 1.6, 5, [25, 60, 5], [25, -60, 5], 180, [0, 0, 1])
