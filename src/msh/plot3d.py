import os
import sys
import math
import numpy as np
import vtk


class Plot3D(object):
    """
    单块结构网格
    """

    def __init__(self, I, J, K, pts):
        self.I, self.J, self.K = I, J, K
        self.X = np.zeros((K, J, I), float)
        self.Y = np.zeros((K, J, I), float)
        self.Z = np.zeros((K, J, I), float)
        self.IBLANK = np.zeros(I * J * K, int)

        for k in range(0, K):
            for j in range(0, J):
                for i in range(0, I):
                    self.X[k][j][i] = pts[k][j][i][0]
                    self.Y[k][j][i] = pts[k][j][i][1]
                    self.Z[k][j][i] = pts[k][j][i][2]

        self.calc_iblank()

    def calc_iblank(self):
        cnt = -1
        for k in range(0, self.K):
            for j in range(0, self.J):
                for i in range(0, self.I):
                    cnt += 1
                    if self.K == 1:
                        if i == 0 or i == self.I - 1 or j == 0 or j == self.J - 1:
                            self.IBLANK[cnt] = 2
                        else:
                            self.IBLANK[cnt] = 1
                    else:
                        if i == self.I - 1 or j == self.J - 1 or k == self.K - 1 or i == 0 or j == 0 or k == 0:
                            self.IBLANK[cnt] = 2
                        else:
                            self.IBLANK[cnt] = 1

    def write_cord(self, cord, fout):
        for k in range(0, self.K):
            for j in range(0, self.J):
                for i in range(0, self.I):
                    if i == 0:
                        fout.write("\n{}".format(cord[k][j][i]))
                    else:
                        fout.write(" {}".format(cord[k][j][i]))

    def write_iblank(self, fout):
        ci = -1
        for k in range(0, self.K):
            for j in range(0, self.J):
                for i in range(0, self.I):
                    ci += 1
                    if i == 0:
                        fout.write("\n{}".format(self.IBLANK[ci]))
                    else:
                        fout.write(" {}".format(self.IBLANK[ci]))

    def output(self, filename):
        fout = open(filename, 'w')

        # BLOCK_NUM
        # I, J, K
        fout.write("{}\n{} {} {}".format(1, self.I, self.J, self.K))

        # X,Y,Z Coordinates
        self.write_cord(self.X, fout)
        self.write_cord(self.Y, fout)
        self.write_cord(self.Z, fout)

        # IBLANK
        self.write_iblank(fout)

        fout.write("\n")
        fout.close()

    @classmethod
    def view(cls, filename):
        reader = vtk.vtkMultiBlockPLOT3DReader()
        reader.BinaryFileOff()
        reader.MultiGridOn()
        reader.SetXYZFileName(filename)
        reader.Update()
        output = reader.GetOutput().GetBlock(0)

        filter = vtk.vtkStructuredGridGeometryFilter()
        filter.SetInputData(output)
        filter.SetExtent(0, 100, 0, 100, 0, 100)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(filter.GetOutputPort())
        # mapper.SetScalarRange(output.GetPointData().GetRange())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        ren = vtk.vtkRenderer()
        ren.AddActor(actor)  # Assign actor to the renderer

        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)  # Create a rendering window and renderer

        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)  # Create a render window interactor
        iren.Initialize()  # Enable user interface interactor

        renWin.Render()
        iren.Start()


R_MIN = 10
R_MAX = 100
U = 61
V = 16
W = 20
fn = "../../result/cylinder.xyz"

if __name__ == "__main__":

    pts = np.zeros((W, V, U, 3), float)
    for w in range(0, W):
        for v in range(0, V):
            for u in range(0, U):
                pts[w][v][u][0] = u
                pts[w][v][u][1] = v
                pts[w][v][u][2] = w

    msh = Plot3D(U, V, W, pts)
    msh.output(fn)

    Plot3D.view(fn)
