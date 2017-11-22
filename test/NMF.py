import unittest
import numpy as np
from src.grid import NMFEntry


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
            self.assertEqual(entry.pri_node_num, ans[k][0])
            self.assertEqual(entry.sec_node_num, ans[k][1])
            self.assertEqual(entry.node_num, ans[k][2])
            self.assertEqual(entry.face_num, ans[k][3])

    def test_real_logic_conversion(self):
        entry = NMFEntry.single('FAR', 10, 4, 5, 33, 4, 17)
        entry.B1Shape = np.array([30, 40, 20, 3])

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

    def test_opposite_logic(self):
        pass


if __name__ == '__main__':
    unittest.main()
