import sys

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage python convert.py foil_file')
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
        f.write("{} {} {}\n".format(cc[0], cc[1], 0.0))
    f.close()
