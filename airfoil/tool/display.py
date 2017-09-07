from matplotlib import pyplot as plt
import sys

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage \'python display.py foil.dat\'")
        exit(-1)

    fn = sys.argv[1]
    pts = []
    fin = open(fn)
    for line in fin:
        x, y, z = line.split()
        pts.append((float(x), float(y), float(z)))
    fin.close()

    (px, py, pz) = zip(*pts)
    plt.plot(px, py)
    plt.gca().set_aspect('equal')
    plt.show()
