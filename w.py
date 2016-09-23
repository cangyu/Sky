class Point_3D(object):
    '3维空间点'

    def __init__(self, _x, _y, _z):
        self.x = _x
        self.y = _y
        self.z = _z


class Airfoil(object):
    '2维翼型描述, 百分比形式'

    def __init__(self, _x, _y_up, _y_down):
        self.x = _x
        self.y_up = _y_up
        self.y_down = _y_down


class Section(object):
    '组成3维机翼的真实剖面'

    def __init__(self, _airfoil, _z, _x_front, _x_tail, _dy, _theta):
        '''
        :param _airfoil: 剖面上的翼型
        :param _z: 剖面与根部的距离，百分比形式
        :param _x_front: 前缘点在x方向上的位置
        :param _x_tail: 后缘点在x方向上的位置
        :param _dy: 由于机翼上反在y方向上引起的偏移
        :param _theta: 绕前缘的扭转角
        '''
        self.airfoil=_airfoil
        self.z=_z
        self.x_front=_x_front
        self.x_tail=_x_tail
        self.dy=_dy
        self.theta=_theta
        self.up_pts=[]
        self.down_pts=[]


