import os
import sys

XFOIL_PATH = os.path.join(os.getcwd(), 'XFOIL6.99', 'xfoil')


def find_alpha(foil, re, ma, cl, iter_cnt=5000):
    foil_path = os.path.join(os.getcwd(), '..', foil + '.dat')

    polar_fn = '_polar.dat'
    polar_path = os.path.join(os.getcwd(), polar_fn)

    '''Generate command file'''
    cmd_fn = '_command.in'
    cmd_path = os.path.join(os.getcwd(), cmd_fn)
    cmd_stream = open(cmd_path, 'w')

    def issue_cmd(cmd):
        cmd_stream.write(cmd + '\n')
        # print(cmd)

    issue_cmd('load ' + foil_path)
    issue_cmd(foil)
    issue_cmd('panel')
    issue_cmd('oper')
    issue_cmd('visc {}'.format(re))
    issue_cmd('M {}'.format(ma))
    issue_cmd('type 1')
    issue_cmd('pacc')
    issue_cmd(polar_fn)
    issue_cmd('')
    issue_cmd('iter')
    issue_cmd(str(iter_cnt))
    issue_cmd('cl {}'.format(cl))
    issue_cmd('')
    issue_cmd('')
    issue_cmd('quit')
    cmd_stream.close()

    '''Execute XFOIL commands'''
    os.system(XFOIL_PATH + ' < ' + cmd_path)
    os.remove(cmd_path)

    '''Extract results'''
    polar_stream = open(polar_path)
    lines = polar_stream.readlines()
    polar_stream.close()
    os.remove(polar_path)

    '''Output'''
    print('\nConfiguration for {} with Re={}, Mach={}, Cl={}:'.format(foil, re, ma, cl))
    for i in [-3, -2, -1]:
        print('\t' + lines[i])


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage \'python find_alpha.py airfoil reynolds mach target_cl\'")
        exit(-1)

    airfoil = sys.argv[1]
    reynolds = int(float(sys.argv[2]))
    mach = float(sys.argv[3])
    target_cl = float(sys.argv[4])

    find_alpha(airfoil, reynolds, mach, target_cl)
