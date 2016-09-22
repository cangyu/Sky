import numpy as np
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot


airfoil_file = "./data/table11.dat"
airfoil_fig = "./fig11.png"
airfoil = open(airfoil_file)

# 平均气动弦长(Mean Aerodynamic Chord)
MAC = 650.0

# 展长
Span = 4270.0
Span_Half = Span / 2

# 3个数组，描述翼型
x = []
y_up = []
y_down = []

for line in airfoil:
    (_x, _y_up, _y_down) = line.split()
    x.append(float(_x) * MAC)
    y_up.append(float(_y_up) * MAC)
    y_down.append(float(_y_down) * MAC)

n = len(x)
x_horizontal_start = []
x_horizontal_end = []
for i in range(0, n):
    x_horizontal_start.append([x[i], (y_up[i]+y_down[i])/2, 0])
    x_horizontal_end.append([x[i],(y_up[i]+y_down[i])/2,Span_Half])

faces = []

# 上表面
section_start_up = []
section_end_up = []
for i in range(0, len(x)):
    section_start_up.append([x[i], y_up[i], 0])
    section_end_up.append([x[i], y_up[i], Span_Half])

for i in range(0, n - 1):
    faces.append([section_start_up[i], section_start_up[i + 1], section_end_up[i]])
    faces.append([section_end_up[i], section_end_up[i + 1], section_start_up[i + 1]])

# 下表面
section_start_down = []
section_end_down = []
for i in range(0, len(x)):
    section_start_down.append([x[i], y_down[i], 0])
    section_end_down.append([x[i], y_down[i], Span_Half])

for i in range(0, n - 1):
    faces.append([section_start_down[i], section_start_down[i + 1], section_end_down[i]])
    faces.append([section_end_down[i], section_end_down[i + 1], section_start_down[i + 1]])

# 前端面
faces.append([x_horizontal_start[0], x_horizontal_start[1], section_start_up[1]])
faces.append([x_horizontal_start[0], x_horizontal_start[1], section_start_down[1]])
for i in range(1, n - 1):
    faces.append([x_horizontal_start[i], section_start_up[i], section_start_up[i + 1]])
    faces.append([x_horizontal_start[i], section_start_down[i], section_start_down[i + 1]])
    faces.append([x_horizontal_start[i], x_horizontal_start[i+1], section_start_up[i + 1]])
    faces.append([x_horizontal_start[i], x_horizontal_start[i+1], section_start_down[i + 1]])
    
# 后端面
faces.append([x_horizontal_end[0], x_horizontal_end[1], section_end_up[1]])
faces.append([x_horizontal_end[0], x_horizontal_end[1], section_end_down[1]])
for i in range(1, n - 1):
    faces.append([x_horizontal_end[i], section_end_up[i], section_end_up[i + 1]])
    faces.append([x_horizontal_end[i], section_end_down[i], section_end_down[i + 1]])
    faces.append([x_horizontal_end[i], x_horizontal_end[i+1], section_end_up[i + 1]])
    faces.append([x_horizontal_end[i], x_horizontal_end[i+1], section_end_down[i + 1]])

#后缘
faces.append([section_start_up[n-1],section_start_down[n-1],section_end_up[n-1]])
faces.append([section_end_up[n-1],section_end_down[n-1],section_start_down[n-1]])

# 简单直机翼
face_num_up = 2 * (n - 1)
face_num_down = 2 * (n - 1)
face_num_start = 4 * (n - 1) - 2
face_num_end = 4 * (n - 1) - 2
face_num_tail = 2

face_num = face_num_up + face_num_down + face_num_start + face_num_end + face_num_tail

wing = mesh.Mesh(np.zeros(face_num, dtype=mesh.Mesh.dtype))
for i in range(0, face_num):
    wing.vectors[i]=faces[i]

wing.save('wing.stl')

#显示
fig=pyplot.figure()
axes=mplot3d.Axes3D(fig)

axes.add_collection3d(mplot3d.art3d.Poly3DCollection(wing.vectors))

scale=wing.points.flatten(-1)
axes.auto_scale_xyz(scale,scale,scale)

#pyplot.show()
