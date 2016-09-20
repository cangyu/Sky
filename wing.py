import matplotlib.pyplot as plt

airfoil_file="./data/table11.dat"
airfoil_fig="./fig11.png"
airfoil=open(airfoil_file)

#3个数组，描述翼型
x=[]
y_up=[]
y_down=[]


for line in airfoil:
    (_x,_y_up,_y_down)=line.split()
    x.append(float(_x))
    y_up.append(float(_y_up))
    y_down.append(float(_y_down))

plt.plot(x, y_up, label="Up")
plt.plot(x, y_down, label="Down")
plt.legend()
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
plt.savefig(airfoil_fig)
plt.close()
airfoil.close()
