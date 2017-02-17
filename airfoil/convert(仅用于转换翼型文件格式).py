import os
import sys
import math

original_filename='naca0012.dat'
converted_filename='naca0012_converted.dat'

airfoil = open(original_filename)

x_up=[]
x_down=[]
y_up=[]
y_down=[]

pts_num=66
cur_line_index=0

for line in airfoil:
    _x, _y = line.split()
    _x=float(_x)
    _y=float(_y)

    if cur_line_index<pts_num:
        x_up.append(_x)
        y_up.append(_y)
    else:
        x_down.append(_x)
        y_down.append(_y)

    cur_line_index+=1

airfoil.close()

for i in range(0, pts_num):
    if x_up[i]!=x_down[i]:
        exit(-1)

airfoil_conv=open(converted_filename,'w')

for i in range(0, pts_num):
    airfoil_conv.write(str(x_up[i])+' '+str(y_up[i])+' '+str(y_down[i])+'\n')

airfoil_conv.close()
