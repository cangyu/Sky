import sys
import numpy as np

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage \'python reverse_order.py foil.dat\'")
        exit(-1)

    fn = sys.argv[1]
    crd = [[], []]
    side = 0

    f = open(fn)
    for k, line in enumerate(f):
        if line.startswith('\n'):
            side = 1
            continue

        x, y = line.split()
        crd[side].append(np.array([float(x), float(y)]))
    f.close()

    f = open(fn, 'w')
    i = len(crd[0]) - 1
    while i >= 0:
        f.write("{:10.6f}\t{:10.6f}\t{:10.6f}\n".format(crd[0][i][0], crd[0][i][1], 0.))
        i -= 1

    for i in range(1, len(crd[1])):
        f.write("{:10.6f}\t{:10.6f}\t{:10.6f}\n".format(crd[1][i][0], crd[1][i][1], 0.))
        
    f.close()
