from ..iges.iges_core import Entity


class Entity116(Entity):
    def __init__(self, _x, _y, _z, _ptr=0):
        """
        Point Entity
        :param _x: X-Coordinate.
        :type _x: float
        :param _y: Y-Coordinate.
        :type _y: float
        :param _z: Z-Coordinate.
        :type _z: float
        :param _ptr: Pointer to the DE of the Sub-figure Definition Entity specifying the display symbol or zero. If zero, no display symbol is specified.
        :type _ptr: int
        """

        super(Entity116, self).__init__(116)
        self.X = float(_x)
        self.Y = float(_y)
        self.Z = float(_z)
        self.PTR = int(_ptr)

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        param = "{},".format(self.directory.entity_type_number)
        param += "{},{},{},{};".format(self.X, self.Y, self.Z, self.PTR)

        return self.to_formatted(param)
