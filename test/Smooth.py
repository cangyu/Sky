import unittest
import numpy as np
import math
from src.grid import LinearTFI2D, uniform, Plot3DBlock, Plot3D, Laplace3D


class EllipticTestCase(unittest.TestCase):
    def test_3d_laplace(self):
        # Delta, R1, R2, U, V, W
        data = [(-10, 4, 25, 16, 41, 21),
                (-1, 2, 5, 16, 33, 16)]

        tfi_grid = []
        laplace_grid = []
        for dt in data:
            tfi = LinearTFI2D(lambda u: np.array([(dt[0] + dt[1]) * (1 - u) + dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[1] * math.cos(math.pi * v) + dt[0], dt[1] * math.sin(math.pi * v), 0]),
                              lambda u: np.array([(dt[0] - dt[1]) * (1 - u) - dt[2] * u, 0, 0]),
                              lambda v: np.array([dt[2] * math.cos(math.pi * v), dt[2] * math.sin(math.pi * v), 0]))
            tfi.calc_grid(uniform(dt[3]), uniform(dt[4]))
            extruded_grid = np.empty((dt[3], dt[4], dt[5], 3))
            for i in range(dt[3]):
                for j in range(dt[4]):
                    for k in range(dt[5]):
                        extruded_grid[i][j][k] = tfi.grid[i][j]
                        extruded_grid[i][j][k][2] = k
            laplace_blk = Laplace3D(extruded_grid)
            laplace_blk.smooth()
            tfi_grid.append(extruded_grid)
            laplace_grid.append(laplace_blk.grid)

        for i in range(len(data)):
            msh = Plot3D()
            msh.add(Plot3DBlock.construct_from_array(tfi_grid[i]))
            msh.save('test_3d_eccentric-{}_tfi.xyz'.format(i))
            msh.clear()
            msh.add(Plot3DBlock.construct_from_array(laplace_grid[i]))
            msh.save('test_3d_eccentric-{}_laplace.xyz'.format(i))
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
