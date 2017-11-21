import unittest
from src.iges import Model, Entity116, Entity110


class IGESTestCase(unittest.TestCase):
    def test_size(self):
        model = Model()
        self.assertTrue(model.size == 0)
        model.add(Entity116(3.14, -2.718, 0.618))
        self.assertTrue(model.size == 1)
        model.add(Entity110((0, 0, 0), (10, 20, 30)))
        self.assertTrue(model.size == 2)

    def test_clear(self):
        model = Model()
        self.assertTrue(model.size == 0)
        model.add(Entity116(3.14, -2.718, 0.618))
        model.add(Entity110((0, 0, 0), (10, 20, 30)))
        self.assertTrue(model.size == 2)
        model.clear()
        self.assertTrue(model.size == 0)

    def test_save(self):
        model = Model()
        model.add(Entity116(3.14, 2.718, 0.618))
        model.save('test_save1.igs')
        self.assertTrue(model.size == 1)
        model.save('test_save2.igs')
        self.assertTrue(model.size == 1)

    def test_pnt(self):
        model = Model()
        pnt = Entity116(3.14, -2.718, 0.618)
        model.add(pnt)
        self.assertTrue(model.size == 1)
        model.save('test_pnt.igs')

    def test_line(self):
        model = Model()
        line = [((0, 0, 0), (10, 20, 30)),
                ((5, 5, 5), (-2, -3.14, 1.618))]
        for l in line:
            model.add(Entity110(l[0], l[1]))
        model.save('test_line.igs')
        self.assertTrue(model.size == len(line))


if __name__ == '__main__':
    unittest.main()
