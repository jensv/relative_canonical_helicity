"""
Created on Feb 16 2017

@author: Jens von der Linden
"""

import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

def get_data(vtk_file):
    r"""
    Get data from a vtk rectilinear grid file.

    Parameters
    ----------
    vtk_file: string
        file name

    Returns
    -------
    data: RectilinearGrid object
        contains all data and grid points.
    """
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_file)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    return data

def get_points(data):
    r"""
    Get points of unstructured grid.
    """
    vtk_points = data.GetPoints()
    number = int(vtk_points.GetNumberOfPoints())
    points = []
    for i in xrange(number):
        points.append(vtk_points.GetPoint(i))
    points = np.asarray(points)
    return points

def get_values(data):
    r"""
    Get values from unstructured grid.
    """
    point_data = data.GetPointData()
    values = vtk_to_numpy(point_data.GetScalars()).astype('float64')
    return values

def read_unstructured_vtk(vtk_file):
    r"""
    Read unstructured grid vtk file
    and return points and values.
    """
    data = get_data(vtk_file)
    points = get_points(data)
    values = get_values(data)
    return points, values
