import os
import sys

NACA456_FILENAME = 'naca456.out'
NACA456_PATH = os.path.join(os.getcwd(), NACA456_FILENAME)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python3 gen.py foil.nml')
        exit(-1)

    foil_desc_path = sys.argv[1]
    foil_name, ext = foil_desc_path.split('.')

    input_path = '_in.dat'
    output_path = '_out.dat'

    f = open(input_path, 'w')
    f.write(foil_desc_path)
    f.close()

    os.system(NACA456_PATH + ' < ' + input_path + ' > ' + output_path)
    os.remove(input_path)
    os.remove(output_path)

    crd = [[], []]
    side = 0
    f = open('naca.gnu')
    for k, line in enumerate(f):
        if len(line.strip()) == 0:
            side = 1
            continue
        x, y = line.split()
        crd[side].append([float(x), float(y)])
    f.close()

    os.remove('naca.gnu')
    os.remove('naca.out')
    os.remove('naca.dbg')

    crd[0].pop(-1)
    crd[1].pop(-1)
    crd[0][-1][0] = 1.0
    crd[1][-1][0] = 1.0

    f = open(foil_name + '.dat', 'w')
    i = len(crd[0]) - 1
    while i >= 0:
        f.write("{:10.6f}\t{:10.6f}\n".format(crd[0][i][0], crd[0][i][1]))
        i -= 1
    for i in range(1, len(crd[1])):
        f.write("{:10.6f}\t{:10.6f}\n".format(crd[1][i][0], crd[1][i][1]))
    f.close()
