import os

if __name__ == '__main__':
    for f in os.listdir(os.getcwd()):
        name, ext = f.split('.')
        if ext in ['nml', 'dat', 'out', 'mod']:
            os.remove(f)
