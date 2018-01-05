import os

rm = ['xyz', 'png', 'dat']

for f in os.listdir('.'):
    fn, ext = f.split('.')
    if ext in rm:
        os.remove(f)
