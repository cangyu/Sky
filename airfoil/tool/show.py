from matplotlib import pyplot as plt
import sys
import os

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage \'python show.py foil.dat\'")
        exit(-1)

    fn = sys.argv[1]
    pts = []
    if fn.endswith('.dat'):
        fin = open(fn)
    else:
        fin = open(os.path.join('..', 'database', fn + '.dat'))
    for line in fin:
        x, y = line.split()
        pts.append((float(x), float(y)))
    fin.close()

    (px, py) = zip(*pts)
    plt.plot(px, py)
    plt.gca().set_aspect('equal')
    plt.show()
