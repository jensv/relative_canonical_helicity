import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def read_vector(vtk_file, component_names):
    r"""
    """
    data = get_data(vtk_file)
    mesh, shape = get_mesh(data)
    vector = get_vector(data, shape, component_names[0],
                        component_names[1], component_names[2])
    return mesh, vector

def read_scalar(vtk_file, component_name):
    r"""
    """
    data = get_data(vtk_file)
    mesh, shape = get_mesh(data)
    scalar = get_scalar(data, shape, component_name)
    return mesh, scalar

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
    reader = vtk.vtkRectilinearGridReader()
    reader.SetFileName(vtk_file)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    data = reader.GetOutput()
    return data


def get_mesh(data):
    r"""
    Reconstruct vtk rectilinear mesh with numpy.

    Parameters
    ----------
    data: RectilinearGrid object
         contains all data and grid points.

    Returns
    -------
    mesh: array of ndarray
        mesh of the rectilinear grid
    shape: tuple
        mesh shape
    """
    x_coord = vtk_to_numpy(data.GetXCoordinates())
    y_coord = vtk_to_numpy(data.GetYCoordinates())
    z_coord = vtk_to_numpy(data.GetZCoordinates())
    mesh = np.meshgrid(x_coord, y_coord, z_coord)
    shape = mesh[0].shape
    return mesh, shape


def get_scalar(data, shape, scalar_name):
    r"""
    Return numpy array of a scalar in vtk RectilinearGrid object.

    Parameters
    ----------
    data: RectilinearGrid object
         contains all data and grid points.
    shape: tuple
        mesh shape
    scalar_name: string
        name of scalar
    Returns
    -------
    scalar: ndarray
        scalar data resized to shape
    """
    shape = (shape[2], shape[0], shape[1])
    data = data.GetPointData()
    scalar = vtk_to_numpy(data.GetArray(scalar_name)).astype('float64')
    scalar.resize(shape)
    scalar = np.swapaxes(np.swapaxes(scalar, 1, 0), 2, 1)
    return scalar


def get_vector(data, shape, x_name,
               y_name, z_name):
    r"""
    Return numpy array of a scalar in vtk RectilinearGrid object.

    Parameters
    ----------
    data: RectilinearGrid object
        contains all data and grid points.
    shape: tuple
        mesh shape
    x_name: string
        name of x component of vector
    y_name: string
        name of y component of vector
    z_name: string
        name of z component of vector
    Returns
    -------
    vector: tuple of ndarray
        vector data resized to shape
    """
    vector_x = get_scalar(data, shape, x_name)
    vector_y = get_scalar(data, shape, y_name)
    vector_z = get_scalar(data, shape, z_name)
    return [vector_x, vector_y, vector_z]
