import subprocess
import os

if __name__ == '__main__':
    p = subprocess.Popen('gfortran ../src/naca456.f90 -o naca456.out', shell=True)
    p.wait()

    for f in os.listdir(os.getcwd()):
        if f.endswith('mod'):
            os.remove(f)
