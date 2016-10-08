import math
import numpy as np

'''剖面形状'''
DESIGN_NUM = 100
RADIUS = 1

x = np.linspace(0, 1, DESIGN_NUM)
y_up = np.linspace(RADIUS, 0, DESIGN_NUM)
y_down = np.linspace(-RADIUS, 0, DESIGN_NUM)

for i in range(0, DESIGN_NUM):
    y_up[i] = math.sqrt(1 - math.pow(x[i], 2))
    y_down[i] = -y_up[i]

fuselage_section = open('./data/fuselage_section01.dat', 'w')
for i in range(0, DESIGN_NUM):
    fuselage_section.write(str(x[i]) + '\t' + str(y_up[i]) + '\t' + str(y_down[i]) + '\n')

fuselage_section.close()

'''机身描述'''
SECTION_NUM=10
FUSELAGE_LEN=6380

frame_shape=['./data/fuselage_section01.dat']*SECTION_NUM
height_dist=[0]*SECTION_NUM
radius_dist=[650]*SECTION_NUM
z_dist=np.linspace(0, FUSELAGE_LEN, SECTION_NUM)

fuselage_desc=open('./data/fuselage1.dat','w')
for i in range(0, SECTION_NUM):
    fuselage_desc.write(frame_shape[i]+'\t'+str(height_dist[i])+'\t'+str(radius_dist[i])+'\t'+str(z_dist[i])+'\n')
fuselage_desc.close()
