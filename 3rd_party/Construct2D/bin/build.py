import os
import shutil

makefile_ref_name = 'Makefile_Linux_MacOSX'
makefile_name = 'Makefile'

if __name__ == '__main__':
    shutil.copyfile(makefile_ref_name, makefile_name)
    os.system('make')
    os.system('make clean')
    os.remove(makefile_name)
