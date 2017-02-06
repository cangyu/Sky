import math
import stl
from stl import mesh
import numpy

def move(_solid, axis, dist):
    if axis == 'x':
        items = [0, 3, 6]
    elif axis == 'y':
        items = [1, 4, 7]
    elif axis == 'z':
        items = [2, 5, 8]
    for p in _solid.points:
        # point items are ((x, y, z), (x, y, z), (x, y, z))
        for i in range(3):
            p[items[i]] += dist


def mirror(_solid, axis):
    if axis == 'x':
        items = [0, 3, 6]
    elif axis == 'y':
        items = [1, 4, 7]
    elif axis == 'z':
        items = [2, 5, 8]
    for p in _solid.points:
        # point items are ((x, y, z), (x, y, z), (x, y, z))
        for i in range(3):
            p[items[i]] = - p[items[i]]

def rotate(_solid, axis, angle):

    c=math.cos(math.radians(angle))
    s=math.sin(math.radians(angle))

    '''旋转矩阵'''
    if axis == 'x':
        M = [[1, 0, 0], [0, c, -s], [0, s, c]]
    elif axis == 'y':
        M = [[c, 0, s], [0, 1, 0], [-s, 0, c]]
    elif axis == 'z':
        M = [[c, -s, 0], [s, c, 0], [0, 0, 1]]

    for p in _solid.points:
        # point items are ((x, y, z), (x, y, z), (x, y, z))
        for i in range(3):
            x=M[0][0]*p[3*i+0]+M[0][1]*p[3*i+1]+M[0][2]*p[3*i+2]
            y=M[1][0]*p[3*i+0]+M[1][1]*p[3*i+1]+M[1][2]*p[3*i+2]
            z=M[2][0]*p[3*i+0]+M[2][1]*p[3*i+1]+M[2][2]*p[3*i+2]

            p[3*i+0]=x
            p[3*i+1]=y
            p[3*i+2]=z

if __name__ == '__main__':

    '''机身'''
    fuselage = mesh.Mesh.from_file('./result/fuselage.stl')

    '''机翼'''
    wing_left = mesh.Mesh.from_file('./result/wing.stl')
    wing_right= mesh.Mesh.from_file('./result/wing.stl')

    mirror(wing_right, 'z')
    move(wing_left, 'x', 2000)
    move(wing_right, 'x', 2000)

    '''垂尾'''
    vertical_stabilizer = mesh.Mesh.from_file('./result/vertical_stabilizer.stl')
    rotate(vertical_stabilizer, 'x', -90)
    move(vertical_stabilizer, 'x', 6000)

    '''平尾'''
    horizontal_stabilizer_left = mesh.Mesh.from_file('./result/horizontal_stabilizer.stl')
    horizontal_stabilizer_right = mesh.Mesh.from_file('./result/horizontal_stabilizer.stl')

    mirror(horizontal_stabilizer_right, 'z')
    move(horizontal_stabilizer_left, 'x', 6000)
    move(horizontal_stabilizer_right, 'x', 6000)

    '''整机'''
    combined = mesh.Mesh(numpy.concatenate([fuselage.data,
                                            wing_left.data,
                                            wing_right.data,
                                            vertical_stabilizer.data,
                                            horizontal_stabilizer_left.data,
                                            horizontal_stabilizer_right.data]))

    combined.save('./result/aircraft.stl')
