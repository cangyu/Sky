from wing import *
from fuselage import *
from PyQt5.QtWidgets import *
import os, time, vtk


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


# 显示生成的stl模型
class STL_Viewer(object):
    def __init__(self):

        self.reader = vtk.vtkSTLReader()
        self.mapper = vtk.vtkPolyDataMapper()
        self.actor = vtk.vtkActor()
        self.ren = vtk.vtkRenderer()
        self.renWin = vtk.vtkRenderWindow()
        self.iren = vtk.vtkRenderWindowInteractor()

    def show(self, _filename):

        self.reader.SetFileName(_filename)

        if vtk.VTK_MAJOR_VERSION <= 5:
            self.mapper.SetInput(self.reader.GetOutput())
        else:
            self.mapper.SetInputConnection(self.reader.GetOutputPort())

        self.actor.SetMapper(self.mapper)
        self.renWin.AddRenderer(self.ren)  # Create a rendering window and renderer
        self.iren.SetRenderWindow(self.renWin)  # Create a render window interactor
        self.ren.AddActor(self.actor)  # Assign actor to the renderer

        self.iren.Initialize()  # Enable user interface interactor
        self.renWin.Render()
        self.iren.Start()


# 整机
class Aircraft(object):
    def __init__(self):
        Wing.update_airfoil_list()

        # component
        self.fuselage = SimpleFuselage('fuselage.stl')
        self.wing = Wing('wing.stl')
        self.vertical_stabilizer = VerticalStabilizer('vertical_stabilizer.stl')
        self.horizontal_stabilizer = HorizontalStabilizer('horizontal_stabilizer.stl')

        # position of each part relevant to the nose of fuselage
        self.delta_Wing = [2000, 0, 0]
        self.delta_VerticalStabilizer = [5300, 0, 0]
        self.delta_HorizontalStabilizer = [5000, 0, 0]

        # widget
        self.widget = QWidget()
        self.comp_tab = QTabWidget()

        # layout
        self.layout = QVBoxLayout()

        # Viewer
        self.viewer = STL_Viewer()

    def init_widget(self):
        self.fuselage.init_widget()
        self.wing.init_widget()
        self.vertical_stabilizer.init_widget()
        self.horizontal_stabilizer.init_widget()

        self.widget.setLayout(self.layout)

        self.layout.addWidget(self.comp_tab)

        self.comp_tab.addTab(self.fuselage.widget, '机身')
        self.comp_tab.addTab(self.wing.widget, '机翼')
        self.comp_tab.addTab(self.vertical_stabilizer.widget, '垂尾')
        self.comp_tab.addTab(self.horizontal_stabilizer.widget, '平尾')

    def set_pos(self, _wing_pos, _vs_pos, _hs_pos):
        self.delta_Wing = _wing_pos
        self.delta_VerticalStabilizer = _vs_pos
        self.delta_HorizontalStabilizer = _hs_pos

    def update_derived_param(self):
        # first get all changed params
        self.fuselage.update_derived_param()
        self.wing.update_derived_param()
        self.vertical_stabilizer.update_derived_param()
        self.horizontal_stabilizer.update_derived_param()

        # calculate related params
        self.fuselage.HeadingDirCapacity = math.pow(self.fuselage.D, 2) * self.fuselage.L / (
        self.wing.S * 2 * self.wing.Span)
        self.fuselage.PitchDirCapacity = math.pow(self.fuselage.D, 2) * self.fuselage.L / (self.wing.S * self.wing.MAC)
        self.horizontal_stabilizer.V_h = self.horizontal_stabilizer.S / self.wing.S * (
        self.horizontal_stabilizer.X_25 - self.wing.X_25) / self.wing.MAC
        self.vertical_stabilizer.V_v = self.vertical_stabilizer.S / self.wing.S * (
        self.vertical_stabilizer.X_25 - self.wing.X_25) / (2 * self.wing.Span)

        # calculate position
        self.delta_Wing[0] = self.wing.X_25 - 0.25 * self.wing.C_root
        self.delta_Wing[1] = self.wing.dY
        self.delta_Wing[2] = self.wing.dZ

        self.delta_VerticalStabilizer[0] = self.vertical_stabilizer.X_25 - 0.25 * self.vertical_stabilizer.C_root
        self.delta_VerticalStabilizer[1] = self.vertical_stabilizer.dZ
        self.delta_VerticalStabilizer[2] = self.vertical_stabilizer.dY

        self.delta_HorizontalStabilizer[0] = self.horizontal_stabilizer.X_25 - 0.25 * self.horizontal_stabilizer.C_root
        self.delta_HorizontalStabilizer[1] = self.vertical_stabilizer.dY
        self.delta_HorizontalStabilizer[2] = self.vertical_stabilizer.dZ

        # update again
        self.fuselage.update_derived_param()
        self.wing.update_derived_param()
        self.vertical_stabilizer.update_derived_param()
        self.horizontal_stabilizer.update_derived_param()

    def generate(self):
        # create a folder with current timestamp to store generated model
        cur_time = time.strftime(r"%Y-%m-%d_%H-%M-%S", time.localtime())
        cur_folder_path = './result/' + cur_time + '/'
        cur_folder_name = cur_time
        os.mkdir(r'%s/result/%s' % (os.getcwd(), cur_folder_name))

        self.fuselage.filename = cur_folder_path + self.fuselage.filename
        self.wing.filename = cur_folder_path + self.wing.filename
        self.vertical_stabilizer.filename = cur_folder_path + self.vertical_stabilizer.filename
        self.horizontal_stabilizer.filename = cur_folder_path + self.horizontal_stabilizer.filename

        # generate component
        self.fuselage.generate()
        self.wing.generate()
        self.vertical_stabilizer.generate()
        self.horizontal_stabilizer.generate()

        # import component
        fuselage = mesh.Mesh.from_file(cur_folder_path + 'fuselage.stl')
        wing_left = mesh.Mesh.from_file(cur_folder_path + 'wing.stl')
        wing_right = mesh.Mesh.from_file(cur_folder_path + 'wing.stl')
        hs_left = mesh.Mesh.from_file(cur_folder_path + 'horizontal_stabilizer.stl')
        hs_right = mesh.Mesh.from_file(cur_folder_path + 'horizontal_stabilizer.stl')
        vs_mid = mesh.Mesh.from_file(cur_folder_path + 'vertical_stabilizer.stl')

        # transform to specified position
        mirror(wing_right, 'z')
        mirror(hs_right, 'z')
        rotate(vs_mid, 'x', -90)

        move(wing_left, 'x', self.delta_Wing[0])
        move(wing_left, 'y', self.delta_Wing[1])
        move(wing_left, 'z', self.delta_Wing[2])

        move(wing_right, 'x', self.delta_Wing[0])
        move(wing_right, 'y', self.delta_Wing[1])
        move(wing_right, 'z', self.delta_Wing[2])

        move(hs_left, 'x', self.delta_HorizontalStabilizer[0])
        move(hs_left, 'y', self.delta_HorizontalStabilizer[1])
        move(hs_left, 'z', self.delta_HorizontalStabilizer[2])

        move(hs_right, 'x', self.delta_HorizontalStabilizer[0])
        move(hs_right, 'y', self.delta_HorizontalStabilizer[1])
        move(hs_right, 'z', self.delta_HorizontalStabilizer[2])

        move(vs_mid, 'x', self.delta_VerticalStabilizer[0])
        move(vs_mid, 'y', self.delta_VerticalStabilizer[1])
        move(vs_mid, 'z', self.delta_VerticalStabilizer[2])

        # 整机
        combined = mesh.Mesh(np.concatenate([fuselage.data,
                                             wing_left.data,
                                             wing_right.data,
                                             vs_mid.data,
                                             hs_left.data,
                                             hs_right.data]))

        combined.save(cur_folder_path + 'aircraft.stl')
        self.viewer.show(cur_folder_path + 'aircraft.stl')
