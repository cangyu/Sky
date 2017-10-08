from ..iges.iges_core import Entity


class Entity110(Entity):
    def __init__(self, _p1, _p2, _form=0):
        """
        Line Entity
        :param _p1: Starting point.
        :param _p2: Ending point.
        :param _form: Form number.
        :type _form: int
        """

        super(Entity110, self).__init__(110)
        self.directory.form_number = _form

        self.X1 = float(_p1[0])
        self.Y1 = float(_p1[1])
        self.Z1 = float(_p1[2])
        self.X2 = float(_p2[0])
        self.Y2 = float(_p2[1])
        self.Z2 = float(_p2[2])

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        ret = "{},".format(self.directory.entity_type_number)
        ret += "{},{},{},".format(self.X1, self.Y1, self.Z1)
        ret += "{},{},{};".format(self.X2, self.Y2, self.Z2)

        return self.to_formatted(ret)
