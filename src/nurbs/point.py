import math

class Point2D(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    @classmethod
    def dist(cls, lhs, rhs):
        ans = math.pow(lhs.x - rhs.x, 2)
        ans += math.pow(lhs.y - rhs.y, 2)
        return math.sqrt(ans)
