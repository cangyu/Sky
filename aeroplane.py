import math
import stl
from stl import mesh
import numpy
import time
import os
from wing import Wing
from fuselage import Fuselage


# 沿坐标轴平移
def move(_solid, axis, dist):
    if axis == 'x':
        items = [0, 3, 6]
    elif axis == 'y':
        items = [1, 4, 7]
    elif axis == 'z':
        items = [2, 5, 8]
    for p in _solid.points:
        # point items are ((x, y, z), (x, y, z), (x, y, z))
        for i in range(3):
            p[items[i]] += dist


# 关于坐标轴镜像
def mirror(_solid, axis):
    if axis == 'x':
        items = [0, 3, 6]
    elif axis == 'y':
        items = [1, 4, 7]
    elif axis == 'z':
        items = [2, 5, 8]
    for p in _solid.points:
        # point items are ((x, y, z), (x, y, z), (x, y, z))
        for i in range(3):
            p[items[i]] = - p[items[i]]


# 绕坐标轴旋转
def rotate(_solid, axis, angle):
    c = math.cos(math.radians(angle))
    s = math.sin(math.radians(angle))

    # 旋转矩阵
    if axis == 'x':
        M = [[1, 0, 0], [0, c, -s], [0, s, c]]
    elif axis == 'y':
        M = [[c, 0, s], [0, 1, 0], [-s, 0, c]]
    elif axis == 'z':
        M = [[c, -s, 0], [s, c, 0], [0, 0, 1]]

    for p in _solid.points:
        # point items are ((x, y, z), (x, y, z), (x, y, z))
        for i in range(3):
            x = M[0][0] * p[3 * i + 0] + M[0][1] * p[3 * i + 1] + M[0][2] * p[3 * i + 2]
            y = M[1][0] * p[3 * i + 0] + M[1][1] * p[3 * i + 1] + M[1][2] * p[3 * i + 2]
            z = M[2][0] * p[3 * i + 0] + M[2][1] * p[3 * i + 1] + M[2][2] * p[3 * i + 2]

            p[3 * i + 0] = x
            p[3 * i + 1] = y
            p[3 * i + 2] = z


class Aircraft(object):
    '''整机'''

    def __init__(self, _js, _jy, _cw, _pw):
        '''
        从参数表生成部件
        :param _js: 机身参数列表
        :param _jy: 机翼参数列表
        :param _cw: 垂尾参数列表
        :param _pw: 平尾参数列表
        '''
        self.jishen_param = _js
        self.jiyi_param = _jy
        self.chuiwei_param = _cw
        self.pingwei_param = _pw

        # 以当前时间创建文件夹，存放结果
        self.cur_folder_name = time.strftime(r"%Y-%m-%d_%H-%M-%S", time.localtime())
        os.mkdir(r'%s/result/%s' % (os.getcwd(), self.cur_folder_name))

        # 创建部件
        f = Fuselage('./result/' + self.cur_folder_name + '/fuselage.stl')
        f.set_parameter(_js[0], _js[1], _js[2], _js[3], _js[4], _js[5], _js[6])
        f.generate()

        w = Wing('./result/' + self.cur_folder_name + '/wing.stl')
        w.generate_wing_linear(_jy[0], _jy[1], _jy[2], _jy[3], _jy[4], _jy[5], _jy[6], _jy[7])

        vs = Wing('./result/' + self.cur_folder_name + '/vertical_stabilizer.stl')
        vs.generate_wing_linear(_cw[0], _cw[1], _cw[2], _cw[3], _cw[4], _cw[5], _cw[6], _cw[7])

        hs = Wing('./result/' + self.cur_folder_name + '/horizontal_stabilizer.stl')
        hs.generate_wing_linear(_pw[0], _pw[1], _pw[2], _pw[3], _pw[4], _pw[5], _pw[6], _pw[7])

    def set_pos(self):
        pass

    def generate(self):
        # 机身
        fuselage = mesh.Mesh.from_file('./result/' + self.cur_folder_name + '/fuselage.stl')

        # 机翼
        wing_left = mesh.Mesh.from_file('./result/' + self.cur_folder_name + '/wing.stl')
        wing_right = mesh.Mesh.from_file('./result/' + self.cur_folder_name + '/wing.stl')

        mirror(wing_right, 'z')
        move(wing_left, 'x', 2000)
        move(wing_right, 'x', 2000)

        # 垂尾
        vertical_stabilizer = mesh.Mesh.from_file('./result/' + self.cur_folder_name + '/vertical_stabilizer.stl')
        rotate(vertical_stabilizer, 'x', -90)
        move(vertical_stabilizer, 'x', 6000)

        # 平尾
        horizontal_stabilizer_left = mesh.Mesh.from_file(
            './result/' + self.cur_folder_name + '/horizontal_stabilizer.stl')
        horizontal_stabilizer_right = mesh.Mesh.from_file(
            './result/' + self.cur_folder_name + '/horizontal_stabilizer.stl')

        mirror(horizontal_stabilizer_right, 'z')
        move(horizontal_stabilizer_left, 'x', 6000)
        move(horizontal_stabilizer_right, 'x', 6000)

        # 整机
        combined = mesh.Mesh(numpy.concatenate([fuselage.data,
                                                wing_left.data,
                                                wing_right.data,
                                                vertical_stabilizer.data,
                                                horizontal_stabilizer_left.data,
                                                horizontal_stabilizer_right.data]))

        combined.save('./result/' + self.cur_folder_name + '/aircraft.stl')
