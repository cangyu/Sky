import matplotlib.pyplot as plt
import numpy as np
import time
from stl import mesh

class Frame(object):
    '机身框的描述'

    def __init__(self, _file):
        '''机身剖面的一半，半径为单位1，高度是相对于半径的比值'''

        self.x = []
        self.y_up = []
        self.y_down = []

        frame = open(_file)
        for line in frame:
            (_x, _y_up, _y_down) = line.split()
            self.x.append(float(_x))
            self.y_up.append(float(_y_up))
            self.y_down.append(float(_y_down))
        frame.close()

    def show(self):
        x_reverse = np.zeros(len(self.x))
        for i in range(0, len(x_reverse)):
            x_reverse[i] = -self.x[i]

        # 右侧
        plt.plot(self.x, self.y_up)
        plt.plot(self.x, self.y_down)
        # 左侧
        plt.plot(x_reverse, self.y_up)
        plt.plot(x_reverse, self.y_down)
        plt.show()

class Section(object):
    '机身真实剖面'

    def __init__(self, _frame, _height, _radius, _z):
        '''
        :param _frame: 剖面基准形状
        :param _height: 基线偏移高度（y方向）
        :param _radius: 剖面半径
        :param _z: 沿机身方向的位置（z方向）
        '''
        self.frame = _frame
        self.height = _height
        self.radius = _radius
        self.z = _z
        self.pts_up = np.zeros((len(_frame.x), 3))
        self.pts_down = np.zeros((len(_frame.x), 3))

        self.calc_section_pts()

    def calc_section_pts(self):
        '''计算围成剖面的所有点的3维坐标'''

        for i in range(0, len(self.frame.x)):
            self.pts_up[i][0] = self.pts_down[i][0] = self.radius * self.frame.x[i]

            self.pts_up[i][1] = self.radius * self.frame.y_up[i]
            self.pts_down[i][1] = self.radius * self.frame.y_down[i]

            self.pts_up[i][2] = self.pts_down[i][2] = self.z

        # y方向上的偏移
        for i in range(0, len(self.frame.x)):
            self.pts_up[i][1] += self.height
            self.pts_down[i][1] += self.height

    def generate_face(self):
        '''该剖面自身围成的曲面，主要用于端部'''

        n = len(self.frame.x)
        mid = np.zeros((n, 3))
        for i in range(0, n):
            mid[i][0] = 0.5 * (self.pts_up[i][0] + self.pts_down[i][0])
            mid[i][1] = 0.5 * (self.pts_up[i][1] + self.pts_down[i][1])
            mid[i][2] = 0.5 * (self.pts_up[i][2] + self.pts_down[i][2])

        face = []
        for i in range(0, n - 1):
            face.append([mid[i], self.pts_up[i], self.pts_up[i + 1]])
            face.append([mid[i], self.pts_down[i], self.pts_down[i + 1]])
            face.append([mid[i], mid[i + 1], self.pts_up[i + 1]])
            face.append([mid[i], mid[i + 1], self.pts_down[i + 1]])

        return face

    def generate_surf(cls, start, end):
        '''相邻两个section之间连成曲面'''

        assert len(start.frame.x) == len(end.frame.x)

        n = len(start.frame.x)
        surf = []

        for i in range(0, n - 1):
            surf.append([start.pts_up[i], start.pts_up[i + 1], end.pts_up[i]])
            surf.append([end.pts_up[i], end.pts_up[i + 1], start.pts_up[i + 1]])
            surf.append([start.pts_down[i], start.pts_down[i + 1], end.pts_down[i]])
            surf.append([end.pts_down[i], end.pts_down[i + 1], start.pts_down[i + 1]])

        surf.append([start.pts_up[n - 1], start.pts_down[n - 1], end.pts_up[n - 1]])
        surf.append([end.pts_up[n - 1], end.pts_down[n - 1], start.pts_down[n - 1]])

        return surf

class Fuselage(object):
    '''由若干剖面组成的机身'''

    def __init__(self, _storage_dst):
        self.filename = _storage_dst
        self.section_dist = []

    def set_parameters(self, _frame, _height, _radius, _z):
        self.frame_dist = _frame
        self.height_dist = _height
        self.radius_dist = _radius
        self.z_dist = _z

    def design(self, _description_file):
        desc=open(_description_file)

        _frame_dist=[]
        _height_dist=[]
        _radius_dist=[]
        _z_dist=[]

        for line in desc:
            (frame_path, height, radius, z)=line.split()
            _frame_dist.append(Frame(frame_path))
            _height_dist.append(float(height))
            _radius_dist.append(float(radius))
            _z_dist.append(float(z))

        self.set_parameters(_frame_dist, _height_dist, _radius_dist, _z_dist)

    def create_sections(self):
        for i in range(0, len(self.z_dist)):
            self.section_dist.append(Section(self.frame_dist[i],self.height_dist[i],self.radius_dist[i],self.z_dist[i]))

    def generate_fuselage_stl(self):
        self.create_sections()

        fuselage_surf=[]
        #表面
        for i in range(0, len(self.section_dist)-1):
            sect=self.section_dist[i].generate_surf(self.section_dist[i],self.section_dist[i+1])
            for triangle in sect:
                fuselage_surf.append(triangle)

        #端部与尾部
        for triangle in self.section_dist[0].generate_face():
            fuselage_surf.append(triangle)
        for triangle in self.section_dist[-1].generate_face():
            fuselage_surf.append(triangle)

        #生成STL格式文件
        fuselage = mesh.Mesh(np.zeros(len(fuselage_surf), dtype=mesh.Mesh.dtype))
        for i in range(0, len(fuselage_surf)):
            fuselage.vectors[i] = fuselage_surf[i]

        fuselage.save(self.filename)

if __name__ == '__main__':

    F1=Fuselage('./result/fuselage_'+str(time.time())+'.stl')
    F1.design('./data/fuselage1.dat')
    F1.generate_fuselage_stl()
