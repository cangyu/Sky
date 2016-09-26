import matplotlib.pyplot as plt
import numpy as np
import math
from stl import mesh

class Airfoil(object):
    '2维翼型描述, 弦长为单位1'

    def __init__(self, _file):
        self.x = []
        self.y_up = []
        self.y_down = []

        airfoil = open(_file)
        for line in airfoil:
            (_x, _y_up, _y_down) = line.split()
            self.x.append(float(_x))
            self.y_up.append(float(_y_up))
            self.y_down.append(float(_y_down))
        airfoil.close()

    def show(self):
        plt.plot(self.x, self.y_up, label="UP")
        plt.plot(self.x, self.y_down, label="DOWN")
        plt.show()


class Section(object):
    '组成3维机翼的真实剖面'

    def __init__(self, _airfoil, _z, _x_front, _x_tail, _dy, _theta):
        '''
        :param _airfoil: 剖面上的翼型
        :param _z: 剖面与根部的距离，百分比形式,[0,1]
        :param _x_front: 前缘点在x方向上的位置
        :param _x_tail: 后缘点在x方向上的位置
        :param _dy: 由于机翼上反角导致的在y方向上引起的偏移
        :param _theta: 绕前缘的扭转角
        '''

        self.airfoil = _airfoil
        self.z = _z
        self.x_front = _x_front
        self.x_tail = _x_tail
        self.dy = _dy
        self.theta = _theta

        self.n = len(_airfoil.x)
        self.chord = np.zeros((self.n, 3))  # 剖面弦线上的离散点
        self.up_pts = np.zeros((self.n, 3))  # 剖面上表面的点
        self.down_pts = np.zeros((self.n, 3))  # 剖面下表面的点

    def calc_discrete_pts(self):
        '''计算描述剖面形状的空间坐标点'''

        # 弦长
        chord_len = self.x_tail - self.x_front

        # 计算正常放置时上、下表面上每点的坐标
        for i in range(0, self.n):
            '''x方向'''
            self.up_pts[i][0] = self.x_front + chord_len * self.airfoil.x[i]
            self.down_pts[i][0] = self.up_pts[i][0]
            '''y方向'''
            self.up_pts[i][1] = chord_len * self.airfoil.y_up[i]
            self.down_pts[i][1] = chord_len * self.airfoil.y_down[i]
            '''z方向'''
            self.up_pts[i][2] = self.z
            self.down_pts[i][2] = self.z

        # 计算由于绕前缘扭转带来的改变
        rotate = complex(math.cos(math.radians(self.theta)),math.sin(math.radians(self.theta)))
        for i in range(1, self.n):
            t = complex(self.up_pts[i][0] - self.up_pts[0][0], self.up_pts[i][1] - self.up_pts[0][1]) * rotate
            self.up_pts[i][0] = self.up_pts[0][0] + t.real
            self.up_pts[i][1] = self.up_pts[0][1] + t.imag
        for i in range(1, self.n):
            t = complex(self.down_pts[i][0] - self.down_pts[0][0], self.down_pts[i][1] - self.down_pts[0][1]) * rotate
            self.down_pts[i][0] = self.down_pts[0][0] + t.real
            self.down_pts[i][1] = self.down_pts[0][1] + t.imag

        # 加上由于上反角导致的y方向上的增量
        for i in range(0, self.n):
            self.up_pts[i][1] += self.dy
            self.down_pts[i][1] += self.dy

        # 计算弦线上每点坐标，3维
        for i in range(0, self.n):
            self.chord[i][0] = 0.5 * (self.up_pts[i][0] + self.down_pts[i][0])
            self.chord[i][1] = 0.5 * (self.up_pts[i][1] + self.down_pts[i][1])
            self.chord[i][2] = 0.5 * (self.up_pts[i][2] + self.down_pts[i][2])

    def generate_face(self):
        face=[]
        face.append([self.up_pts[0],self.up_pts[1],self.chord[1]])
        face.append([self.down_pts[0], self.down_pts[1], self.chord[1]])
        for i in range(1, self.n-1):
            face.append([self.chord[i], self.up_pts[i], self.up_pts[i + 1]])
            face.append([self.chord[i], self.down_pts[i], self.down_pts[i + 1]])
            face.append([self.chord[i], self.chord[i + 1], self.up_pts[i + 1]])
            face.append([self.chord[i], self.chord[i + 1], self.down_pts[i + 1]])

        return face

    def generate_surf(cls, start, end):
        '''
        生成两个剖面之间的曲面，两个剖面要有相同的样本点n
        :param start: 起始剖面
        :param end: 结束剖面
        :return: 一个数组，包含了构成曲面的所有三角形
        '''
        assert start.n==end.n

        surf=[]
        n=start.n
        #上表面
        for i in range(0, n - 1):
            surf.append([start.up_pts[i], start.up_pts[i + 1], end.up_pts[i]])
            surf.append([end.up_pts[i], end.up_pts[i + 1], start.up_pts[i + 1]])

        #下表面
        for i in range(0, n - 1):
            surf.append([start.down_pts[i], start.down_pts[i + 1], end.down_pts[i]])
            surf.append([end.down_pts[i], end.down_pts[i + 1], start.down_pts[i + 1]])

        #后缘垂面
        surf.append([start.up_pts[n - 1], start.down_pts[n - 1], end.up_pts[n - 1]])
        surf.append([end.up_pts[n - 1], end.down_pts[n - 1], start.down_pts[n - 1]])

        return surf

class Wing(object):
    '''简单机翼，由一个个剖面组成'''

    def __init__(self, _len):
        self.len=_len

        self.airfoil_list=[]
        self.z_list = []
        self.x_front_list = []
        self.x_tail_list = []
        self.dy_list = []
        self.theta_list = []
        self.section_list=[]

    def calc_sections(self):
        for i in range(0, len(self.z_list)):
            self.section_list.append(Section(self.airfoil_list[i],self.z_list[i]*self.len,self.x_front_list[i],self.x_tail_list[i],self.dy_list[i], self.theta_list[i]))
        for i in range(0, len(self.z_list)):
            self.section_list[i].calc_discrete_pts()


    def generate_wing_stl(self):
        wing_surf=[]

        #上下表面与尾缘
        for i in range(0, len(self.section_list)-1):
            for face in self.section_list[i].generate_surf(self.section_list[i],self.section_list[i+1]):
                wing_surf.append(face)

        #根部与翼梢
        for face in self.section_list[0].generate_face():
            wing_surf.append(face)
        for face in self.section_list[-1].generate_face():
            wing_surf.append(face)

        #生成STL格式文件
        wing = mesh.Mesh(np.zeros(len(wing_surf), dtype=mesh.Mesh.dtype))
        for i in range(0, len(wing_surf)):
            wing.vectors[i] = wing_surf[i]

        wing.save('wing.stl')


if __name__ == '__main__':

    '''
    简单后掠机翼示例
    30度后掠，3度上反，2度负扭转
    '''
    a1 = Airfoil("./data/table11.dat")
    Span=4270/2
    section_num=12
    max_chord_len=650
    min_chord_len=210
    sweep_back_angle=30
    dihedral_angle=0
    twist_angle=0

    w1=Wing(Span)
    #剖面翼型
    for i in range(0, section_num):
        w1.airfoil_list.append(a1)
    #剖面分布
    w1.z_list=np.linspace(0, 1, section_num)
    #前缘点
    w1.x_front_list=np.linspace(0, Span*math.tan(math.radians(sweep_back_angle)), section_num)
    #后缘点
    l=np.linspace(max_chord_len,min_chord_len, section_num)
    for i in range(0, section_num):
        w1.x_tail_list.append(w1.x_front_list[i]+l[i])
    #上反
    w1.dy_list=np.linspace(0, Span*math.tan(math.radians(dihedral_angle)), section_num)
    #扭转
    for i in range(0, section_num):
        w1.theta_list.append(twist_angle)

    w1.calc_sections()
    w1.generate_wing_stl()
