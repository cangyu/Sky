import sys

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage \'python add_zero.py foil.dat\'")
        exit(-1)

    fn = sys.argv[1]
    f = open(fn)
    crd = []
    for line in f:
        (x, y) = line.split()
        crd.append((float(x), float(y)))
    f.close()

    f = open(fn, 'w')
    for cc in crd:
        f.write("{:10.6f}\t{:10.6f}\t{:10.6f}\n".format(cc[0], cc[1], 0.))
    f.close()
