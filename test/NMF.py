import unittest
import numpy as np
import math
from src.grid import NeutralMapFile, NMFEntry, uniform, LinearTFI2D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class NMFTestCase(unittest.TestCase):
    def test_construction(self):
        data = [NMFEntry.single('WALL', 1, 1, 1, 47, 1, 26),
                NMFEntry.single('WALL', 1, 2, 1, 47, 1, 26),
                NMFEntry.one2one(1, 3, 1, 26, 1, 33, 2, 4, 1, 26, 1, 33, False),
                NMFEntry.single('Outflow', 1, 4, 1, 26, 1, 33),
                NMFEntry.one2one(1, 5, 1, 33, 1, 47, 4, 6, 1, 33, 1, 47, False),
                NMFEntry.single('WALL', 1, 6, 1, 33, 1, 47),
                NMFEntry.single('WALL', 2, 1, 1, 19, 1, 26),
                NMFEntry.single('WALL', 2, 2, 1, 19, 1, 26),
                NMFEntry.single('Inflow', 2, 3, 1, 26, 1, 33),
                NMFEntry.one2one(2, 5, 1, 33, 1, 19, 3, 6, 1, 33, 1, 19, False),
                NMFEntry.single('WALL', 2, 6, 1, 33, 1, 19),
                NMFEntry.single('WALL', 3, 1, 1, 19, 1, 24),
                NMFEntry.single('WALL', 3, 2, 1, 19, 1, 24),
                NMFEntry.single('Inflow', 3, 3, 1, 24, 1, 33),
                NMFEntry.one2one(3, 4, 1, 24, 1, 33, 4, 3, 1, 24, 1, 33, False),
                NMFEntry.single('WALL', 3, 5, 1, 33, 1, 19),
                NMFEntry.single('WALL', 4, 1, 1, 47, 1, 24),
                NMFEntry.single('WALL', 4, 2, 1, 47, 1, 24),
                NMFEntry.single('Outflow', 4, 4, 1, 24, 1, 33),
                NMFEntry.single('WALL', 4, 5, 1, 33, 1, 47)]

        # pri_node_num, sec_node_num, node_num, face_num
        ans = np.array([[47, 26, 47 * 26, 46 * 25],
                        [47, 26, 47 * 26, 46 * 25],
                        [26, 33, 26 * 33, 25 * 32],
                        [26, 33, 26 * 33, 25 * 32],
                        [33, 47, 33 * 47, 32 * 46],
                        [33, 47, 33 * 47, 32 * 46],
                        [19, 26, 19 * 26, 18 * 25],
                        [19, 26, 19 * 26, 18 * 25],
                        [26, 33, 26 * 33, 25 * 32],
                        [33, 19, 33 * 19, 32 * 18],
                        [33, 19, 33 * 19, 32 * 18],
                        [19, 24, 19 * 24, 18 * 23],
                        [19, 24, 19 * 24, 18 * 23],
                        [24, 33, 24 * 33, 23 * 32],
                        [24, 33, 24 * 33, 23 * 32],
                        [33, 19, 33 * 19, 32 * 18],
                        [47, 24, 47 * 24, 46 * 23],
                        [47, 24, 47 * 24, 46 * 23],
                        [24, 33, 24 * 33, 23 * 32],
                        [33, 47, 33 * 47, 32 * 46]])

        print('\n' + NMFEntry.NMF_LITERAL)
        for k, entry in enumerate(data):
            print(entry)
            self.assertEqual(entry.pri_node_num(), ans[k][0])
            self.assertEqual(entry.sec_node_num(), ans[k][1])
            self.assertEqual(entry.node_num, ans[k][2])
            self.assertEqual(entry.face_num, ans[k][3])

    def test_real_logic_conversion(self):
        entry = NMFEntry.single('FAR', 10, 4, 5, 33, 4, 17)
        entry.shape_of_blk1 = (30, 40, 20, 3)

        data = [[(0, 0), (29, 4, 3)],
                [(28, 0), (29, 32, 3)],
                [(28, 13), (29, 32, 16)],
                [(0, 13), (29, 4, 16)],
                [(14, 7), (29, 18, 10)]]

        for k, dt in enumerate(data):
            lp = dt[0]
            rp = dt[1]
            self.assertTupleEqual(lp, tuple(entry.real2logic(rp)))
            self.assertTupleEqual(rp, tuple(entry.logic2real(lp)))

    def test_counterpart(self):
        length = 100
        height1 = 40
        height2 = 60
        height3 = 10
        u_num = 6
        v_num = 5

        tfi = LinearTFI2D(lambda u: np.array([u * length, 4 * height3 * u * (1 - u), 0]),
                          lambda v: np.array([0, v * height1, 0]),
                          lambda u: np.array([u * length, (height1 * (1 - u * u) + height2 * u * u), 0]),
                          lambda v: np.array([length, v * height2, 0]))
        tfi.calc_grid(uniform(u_num), uniform(v_num))
        grid = tfi.grid

        plt.scatter(grid[:, :, 0].flatten(), grid[:, :, 1].flatten())
        plt.plot()
        plt.show()

        blk1 = np.empty((12, 6, 5, 3), float)
        blk2 = np.empty((6, 4, 5, 3), float)

        for i in range(12):
            for j in range(6):
                for k in range(5):
                    blk1[i][j][k][0] = math.pi * i
                    blk1[i][j][k][1] = grid[j][k][0]
                    blk1[i][j][k][2] = grid[j][k][1]

        for i in range(6):
            for j in range(4):
                for k in range(5):
                    blk2[i][j][k][0] = math.pi * (16 - i)
                    blk2[i][j][k][1] = grid[4 - k][4 - j][0]
                    blk2[i][j][k][2] = grid[4 - k][4 - j][1]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        xx = blk1[:, :, :, 0].flatten()
        yy = blk1[:, :, :, 1].flatten()
        zz = blk1[:, :, :, 2].flatten()
        ax.scatter(xx, yy, zz, c='b')

        xx = blk2[:, :, :, 0].flatten()
        yy = blk2[:, :, :, 1].flatten()
        zz = blk2[:, :, :, 2].flatten()
        ax.scatter(xx, yy, zz, c='r', marker='^')
        plt.show()

        blk = [blk1, blk2]
        nmf = NeutralMapFile(blk)
        entry = NMFEntry.one2one(1, 4, 1, 5, 2, 5, 2, 4, 1, 4, 1, 5)
        entry.shape_of_blk1 = blk1.shape
        entry.shape_of_blk2 = blk2.shape
        entry.NMF = nmf
        entry.direction_auto_mapping()

        data1 = [[(0, 3), (1, 3), (2, 3), (3, 3), (4, 3)],
                 [(0, 2), (1, 2), (2, 2), (3, 2), (4, 2)],
                 [(0, 1), (1, 1), (2, 1), (3, 1), (4, 1)],
                 [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)]]

        data2 = [[(0, 4), (0, 3), (0, 2), (0, 1), (0, 0)],
                 [(1, 4), (1, 3), (1, 2), (1, 1), (1, 0)],
                 [(2, 4), (2, 3), (2, 2), (2, 1), (2, 0)],
                 [(3, 4), (3, 3), (3, 2), (3, 1), (3, 0)]]

        for i in range(4):
            for j in range(5):
                lp1 = data1[i][j]
                lp2 = data2[i][j]
                self.assertTupleEqual(tuple(entry.counterpart(lp1, 1)), lp2)
                self.assertTupleEqual(tuple(entry.counterpart(lp2, 2)), lp1)


if __name__ == '__main__':
    unittest.main()
