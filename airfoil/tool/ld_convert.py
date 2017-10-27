import sys
from math import sin, cos, radians

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage \'python ld_convert.py r theta\'")
        exit(-1)

    r = float(sys.argv[1])
    theta = radians(float(sys.argv[2]))
    rr = (r*cos(theta)-sin(theta))/(r*sin(theta)+cos(theta))
    print(rr)