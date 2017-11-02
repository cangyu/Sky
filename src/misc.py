import os
import sys
import win32com.client
import subprocess
import math
import numpy as np
from numpy.linalg import norm


def array_smart_copy(src, dst):
    """
    Copy data without overflow.
    :param src: Source.
    :param dst: Destination.
    :return: None.
    """

    n = min(len(src), len(dst))
    for i in range(n):
        dst[i] = src[i]


def normalize(vec):
    """
    Normalize the input vector.
    :param vec: Original vector.
    :return: Normalized vector.
    """

    v = np.copy(vec)
    tmp = norm(v, 2)
    return v if math.isclose(tmp, 0) else v / tmp


def pnt_dist(lhs, rhs):
    """
    Calculate the distance between two points with the great common degree.
    :param lhs: The 1st point.
    :param rhs: The 2nd point.
    :return: The distance between the 2 points.
    :rtype: float
    """

    ans = 0.
    dim = min(len(lhs), len(rhs))

    for i in range(dim):
        ans += math.pow(lhs[i] - rhs[i], 2)

    return math.sqrt(ans)


def view(file):
    app = win32com.client.Dispatch('catia.application')
    app.Visible = True
    app.DisplayFileAlerts = True
    doc = app.Documents.Open(os.path.abspath(file))


def convert(file, target_format):
    cur_dir, filename = os.path.split(file)
    base, ext = os.path.splitext(filename)
    fileout = os.path.join(cur_dir, base + '.' + target_format)

    app = win32com.client.Dispatch('catia.application')
    app.DisplayFileAlerts = False
    doc = app.Documents.Open(os.path.abspath(file))

    doc.ExportData(fileout, format)
    doc.Close()
    app.Quit()
    print("Conversion Done!")

    return fileout


class XFoil(object):
    def __init__(self):
        pass
