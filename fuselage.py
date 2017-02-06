import math
import numpy as np
import matplotlib.pyplot as plt
from stl import mesh

L=6380
L0=0.25*L
L2=0.3*L
L1=L-L0-L2

D=L1/8
h0=0.6*D
h1=D-h0

h2=0.28*D
h3=D-h2

dh=h0-h2

def get_y_up(x):
    if x<0:
        return 0
    elif x<L0:
        return h0*math.sqrt(1-math.pow((x-L0)/L0,2))
    elif x<L0+L1:
        return h0
    elif x<L:
        return dh+h2*math.sqrt(1-math.pow((x-(L0+L1))/L2,2))
    else:
        return dh

def get_y_down(x):
    if x<0:
        return 0
    elif x<L0:
        return -h1*math.sqrt(1-math.pow((x-L0)/L0,2))
    elif x<L0+L1:
        return -h1
    elif x<L:
        return dh-h3*math.sqrt(1-math.pow((x-(L0+L1))/L2,2))
    else:
        return dh


delta_step=5
x_tou=np.arange(0, L0, delta_step)
x_shen=np.arange(L0, L0+L1, delta_step)
x_wei=np.arange(L0+L1, L, delta_step)

x=[]
y_up=[]
y_down=[]

for i in range(0, len(x_tou)):
    x.append(x_tou[i])

for i in range(0, len(x_shen)):
    x.append(x_shen[i])

for i in range(0, len(x_wei)):
    x.append(x_wei[i])

x.append(L)

for i in range(0, len(x)):
    y_up.append(get_y_up(x[i]))
    y_down.append(get_y_down(x[i]))

y_mid=[]

for i in range(0, len(x)):
    y_mid.append((y_up[i]+y_down[i])/2)


plt.plot(x, y_up, label="UP")
plt.plot(x, y_down, label="DOWN")
plt.plot(x, y_mid, label="MID")
plt.gca().set_aspect(1)


fuselage_surf=[]

for i in range(0, len(x)-1):
    rl=(y_up[i]-y_down[i])/2
    rr=(y_up[i+1]-y_down[i+1])/2

    theta=0
    delta_theta=10
    while theta<360:
        pl0=[x[i], rl*math.cos(math.radians(theta)), rl*math.sin(math.radians(theta))]
        pl1=[x[i], rl*math.cos(math.radians(theta+delta_theta)), rl*math.sin(math.radians(theta+delta_theta))]
        pr0 = [x[i+1], rr * math.cos(math.radians(theta)), rr * math.sin(math.radians(theta))]
        pr1 = [x[i+1], rr * math.cos(math.radians(theta + delta_theta)), rr * math.sin(math.radians(theta + delta_theta))]

        fuselage_surf.append([pl0, pl1, pr0])
        fuselage_surf.append([pr0, pr1, pl1])

        theta=theta+delta_theta

fuselage = mesh.Mesh(np.zeros(len(fuselage_surf), dtype=mesh.Mesh.dtype))
for i in range(0, len(fuselage_surf)):
    fuselage.vectors[i] = fuselage_surf[i]

fuselage.save('./result/fuselage.stl')

#plt.show()
