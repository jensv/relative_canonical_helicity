# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 13:48:25 2015

@author: Jens von der Linden
"""

import numpy as np
from pyvisfile.vtk import (write_structured_grid,
                           UnstructuredGrid,
                           DataArray,
                           AppendedDataXMLGenerator,
                           VTK_VERTEX, VF_LIST_OF_VECTORS,
                           VF_LIST_OF_COMPONENTS)
from pytools.obj_array import make_obj_array
from os.path import exists
import scipy.io.idl as idl
from scipy.interpolate import griddata


def plane_points_3d(data, z_position):
    r"""
    Returns 3d points from dictionary of a 2D scan of RSX.
    """
    assert ('x_out' in data.keys() and
            'y_out' in data.keys()), 'No x_out and y_out keys in data '
    points_3d = np.dstack((data['x_out'],
                           data['y_out'],
                           np.ones(data['x_out'].size)*z_position))[0]
    return points_3d


def write_scalar_data_to_vtk(file_name, time_point, z_position, labels,
                             data_dicts):
    r"""
    Writes scalars to an unstrcutured grid VTK file.
    """
    data_points = plane_points_3d(data_dicts[0], z_position)
    if len(data_dicts) > 1:
        msg = 'scalars not measured at same points'
        for data_dict in data_dicts:
            assert np.allclose(data_points,
                               plane_points_3d(data_dict, z_position)), msg
    data = [(labels[i],
             data_dict['a_out'][time_point]) for i, data_dict in enumerate(data_dicts)]
    write_data_to_unstructured_vtk(file_name, data, data_points)


def write_data_to_unstructured_vtk(file_name, data, points):
    r"""
    Writes data to an unstructured grid VTK file.
    Data is a list of points and values, the values can be scalar, or 3d vector.
    """
    n_points = points.shape[0]
    grid = UnstructuredGrid((n_points, DataArray("points",
                                                 points,
                                                 vector_format=VF_LIST_OF_VECTORS)),
                            cells=np.arange(n_points, dtype=np.uint32),
                            cell_types=np.asarray([VTK_VERTEX] * n_points,
                                                  dtype=np.uint8))
    for name, field in data:
        print field.astype('float64')
        grid.add_pointdata(DataArray(name, field.astype('float64'),
                           vector_format=VF_LIST_OF_COMPONENTS))
    if exists(file_name):
        raise RuntimeError("output file '%s' already exists" % file_name)

    outf = open(file_name, "w")
    compressor = None
    AppendedDataXMLGenerator(compressor)(grid).write(outf)
    outf.close()


def write_vector_data_to_vtk(file_name, time_point, z_position, labels,
                             data_dicts):
    r"""
    Writes a vector to an unstructured grid VTK file.
    """
    data_points_x = plane_points_3d(data_dicts[0], z_position)
    data_points_y = plane_points_3d(data_dicts[1], z_position)
    data_points_z = plane_points_3d(data_dicts[2], z_position)
    x_min = np.nanmin(np.concatenate(data_points_x[:][0],
                                     data_points_y[:][0],
                                     data_points_z[:][0]))
    x_max = np.nanmax(np.concatenate(data_points_x[:][0],
                                     data_points_y[:][0],
                                     data_points_z[:][0]))
    y_min = np.nanmin(np.concatenate(data_points_x[:][1],
                                     data_points_y[:][1],
                                     data_points_z[:][1]))
    y_max = np.nanmax(np.concatenate(data_points_x[:][1],
                                     data_points_y[:][1],
                                     data_points_z[:][1]))
    x_points = np.linspace(x_min, x_max, 100)
    y_points = np.linspace(y_min, y_max, 100)
    data_points = np.dstack(x_points, y_points)
    gridx, gridy = np.meshgrid(x_points, y_points)
    interpolated_x = griddata(np.dstack((data_dicts[0]['x_out'],
                                         data_dicts[0]['y_out']))[0],
                              data_dicts[0]['a_out'][time_point],
                              (gridx, gridy))
    interpolated_y = griddata(np.dstack((data_dicts[1]['x_out'],
                                         data_dicts[1]['y_out']))[0],
                              data_dicts[1]['a_out'][time_point],
                              (gridx, gridy))
    interpolated_z = griddata(np.dstack((data_dicts[2]['x_out'],
                                         data_dicts[2]['y_out']))[0],
                              data_dicts[2]['a_out'][time_point],
                              (gridx, gridy))
    interpolated_field = np.dstack((interpolated_x,
                                    interpolated_y,
                                    interpolated_z))[0]
    data = [(labels[0],
             interpolated_field)]
    write_data_to_unstructured_vtk(file_name, data, data_points)
