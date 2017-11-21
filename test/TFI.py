import unittest
import math
import numpy as np
from src.grid import LinearTFI2D, LinearTFI3D, Plot3D, Plot3DBlock
from src.misc import share


class TFITestCase(unittest.TestCase):
    def test_2d_rect(self):
        # L, W, U, V
        data = [(5, 4, 11, 9),
                (8, 8, 31, 21)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([dt[0] * u, 0, 0]),  # C1
                              lambda v: np.array([0, dt[1] * v, 0]),  # C3
                              lambda u: np.array([dt[0] * u, dt[1], 0]),  # C2
                              lambda v: np.array([dt[0], dt[1] * v, 0]))  # C4
            tfi.calc_grid(np.linspace(0, 1, dt[2]),
                          np.linspace(0, 1, dt[3]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_rect-{}.xyz'.format(k))

    def test_2d_circle(self):
        # R1, R2, U, V
        data = [(1, 2, 6, 11),
                (0, 5, 16, 33)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([(1 - u) * dt[0] + u * dt[1], 0, 0]),
                              lambda v: np.array([dt[0] * math.cos(0.5 * math.pi * v), dt[0] * math.sin(0.5 * math.pi * v), 0]),
                              lambda u: np.array([0, (1 - u) * dt[0] + u * dt[1], 0]),
                              lambda v: np.array([dt[1] * math.cos(0.5 * math.pi * v), dt[1] * math.sin(0.5 * math.pi * v), 0]))

            tfi.calc_grid(np.linspace(0, 1, dt[2]),
                          np.linspace(0, 1, dt[3]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_circle-{}.xyz'.format(k))

    def test_2d_eccentric(self):
        # Delta, R1, R2, U, V
        data = [(-10, 4, 25, 16, 41),
                (-10, 0, 25, 16, 41),
                (-1, 2, 5, 16, 33)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([(dt[0] + dt[1]) * (1 - u) + dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[1] * math.cos(math.pi * v) + dt[0], dt[1] * math.sin(math.pi * v), 0]),
                              lambda u: np.array([(dt[0] - dt[1]) * (1 - u) - dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[2] * math.cos(math.pi * v), dt[2] * math.sin(math.pi * v), 0]))

            tfi.calc_grid(np.linspace(0, 1, dt[3]),
                          np.linspace(0, 1, dt[4]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_eccentric-{}.xyz'.format(k))

    def test_2d_crv_rect(self):
        # L, H1, H2, H3
        data = [(100, 40, 60, 10, 50, 25)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI2D(lambda u: np.array([u * dt[0], 4 * dt[3] * u * (1 - u), 0]),
                              lambda v: np.array([0, v * dt[1], 0]),
                              lambda u: np.array([u * dt[0], (dt[1] * (1 - u * u) + dt[2] * u * u), 0]),
                              lambda v: np.array([dt[0], v * dt[2], 0]))

            tfi.calc_grid(np.linspace(0, 1, dt[4]),
                          np.linspace(0, 1, dt[5]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_2d_crv_rect-{}.xyz'.format(k))

    def test_3d_cuboid(self):
        # L, W, H, U, V, W
        data = [(5, 4, 3, 11, 9, 5),
                (10, 10, 10, 21, 31, 41)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI3D(lambda v, w: np.array([0, v * dt[1], w * dt[2]]),
                              lambda v, w: np.array([dt[0], v * dt[1], w * dt[2]]),
                              lambda w, u: np.array([u * dt[0], 0, w * dt[2]]),
                              lambda w, u: np.array([u * dt[0], dt[1], w * dt[2]]),
                              lambda u, v: np.array([u * dt[0], v * dt[1], 0]),
                              lambda u, v: np.array([u * dt[0], v * dt[1], dt[2]]))
            tfi.calc_grid(np.linspace(0, 1.0, dt[3]),
                          np.linspace(0, 1.0, dt[4]),
                          np.linspace(0, 1.0, dt[5]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_3d_cuboid-{}.xyz'.format(k))

    def test_3d_sect(self):
        # R_MIN, R_MAX, THETA_MIN, THETA_MAX, H_MIN, H_MAX, U, V, W
        data = [(5, 20, 60, 120, -50, 50, 31, 16, 61)]
        ans = []

        for k, dt in enumerate(data):
            tfi = LinearTFI3D(lambda v, w: np.array([dt[0] * math.cos(math.radians(share(v, dt[2], dt[3]))), dt[0] * math.sin(math.radians(share(v, dt[2], dt[3]))), share(w, dt[4], dt[5])]),
                              lambda v, w: np.array([dt[1] * math.cos(math.radians(share(v, dt[2], dt[3]))), dt[1] * math.sin(math.radians(share(v, dt[2], dt[3]))), share(w, dt[4], dt[5])]),
                              lambda w, u: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(dt[2])), share(u, dt[0], dt[1]) * math.sin(math.radians(dt[2])), share(w, dt[4], dt[5])]),
                              lambda w, u: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(dt[3])), share(u, dt[0], dt[1]) * math.sin(math.radians(dt[3])), share(w, dt[4], dt[5])]),
                              lambda u, v: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(share(v, dt[2], dt[3]))), share(u, dt[0], dt[1]) * math.sin(math.radians(share(v, dt[2], dt[3]))), dt[4]]),
                              lambda u, v: np.array([share(u, dt[0], dt[1]) * math.cos(math.radians(share(v, dt[2], dt[3]))), share(u, dt[0], dt[1]) * math.sin(math.radians(share(v, dt[2], dt[3]))), dt[5]]))
            tfi.calc_grid(np.linspace(0, 1.0, dt[6]),
                          np.linspace(0, 1.0, dt[7]),
                          np.linspace(0, 1.0, dt[8]))
            ans.append(tfi.grid)

        self.assertTrue(len(data) == len(ans))

        msh = Plot3D()
        for k, g in enumerate(ans):
            msh.clear()
            blk = Plot3DBlock.construct_from_array(g)
            msh.add(blk)
            msh.save('test_3d_sect-{}.xyz'.format(k))


if __name__ == '__main__':
    unittest.main()
