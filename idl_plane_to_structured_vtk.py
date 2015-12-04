# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 11:14:38 2015

@author: Jens von der Linden
"""

import numpy as np
from pyvisfile.vtk import write_structured_grid
from pytools.obj_array import make_obj_array
from os.path import exists
import scipy.io.idl as idl
from scipy.interpolate import griddata

def read_data(file_name):
    r"""
    Read IDL sav file.
    """
    data_dict = idl.readsav(file_name)
    return data_dict


def resample_on_structutred_grid(data_dict,
                                 time_point,
                                 x_min, x_max,
                                 y_min, y_max,
                                 x_points=50,
                                 y_points=50,
                                 to_clip=[0, -1, 0, -1],
                                 method='linear'):
    r"""
    Resample quantity from measurement grid to structured grid.
    """
    x_points = np.linspace(x_min, x_max, 100)
    y_points = np.linspace(y_min, y_max, 100)
    x_grid, y_grid = np.meshgrid(x_points, y_points)
    quantity_interpolated = griddata(np.dstack((data_dict['x_out'],
                                                data_dict['y_out']))[0],
                                                data_dict['a_out'][time_point],
                                                (x_grid, y_grid),
                                                method=method)
    x_slice = slice(to_clip[0], to_clip[1])
    y_slice = slice(to_clip[2], to_clip[3])
    quantity_interpolated = quantity_interpolated[x_slice, y_slice]
    x_grid = x_grid[x_slice, y_slice]
    y_grid = y_grid[x_slice, y_slice]
    return quantity_interpolated, x_grid, y_grid


def determine_sample_bounds(data_dicts):
    r"""
    Determine the x,y bounds of the 2D sampling space.
    """
    x_mins = np.asarray([np.nanmin(data_dict['x_out']) for data_dict in data_dicts])
    x_maxs = np.asarray([np.nanmax(data_dict['x_out']) for data_dict in data_dicts])
    y_mins = np.asarray([np.nanmin(data_dict['y_out']) for data_dict in data_dicts])
    y_maxs = np.asarray([np.nanmax(data_dict['y_out']) for data_dict in data_dicts])
    x_min = x_mins.max()
    x_max = x_maxs.min()
    y_min = y_mins.max()
    y_max = y_maxs.min()
    return x_min, x_max, y_min, y_max


def resample_scalar(scalar_dict,
                    time_point,
                    x_min, x_max,
                    y_min, y_max,
                    **kwargs):
    r"""
    Resample scalar from idl readout from measurement grid to same structured
    grid.
    """
    x_min, x_max, y_min, y_max = determine_sample_bounds([scalar_dict])
    (quanitity_interpolated,
     x_grid, y_grid) = resample_on_structutred_grid(scalar_dict,
                                                    time_point,
                                                    x_min, x_max,
                                                    y_min, y_max,
                                                    **kwargs)
    return quanitity_interpolated, x_grid, y_grid


def resample_vector(vector_x_dict,
                    vector_y_dict,
                    vector_z_dict,
                    time_point,
                    x_min, x_max,
                    y_min, y_max,
                    **kwargs):
    r"""
    Resample vector from idl readout from measurement grid to same structured
    grid.
    """
    data_dicts = [vector_x_dict, vector_y_dict, vector_z_dict]
    x_min, x_max, y_min, y_max = determine_sample_bounds(data_dicts)
    vector_interpolated = []
    for data_dict in data_dicts:
        (quanitity_interpolated,
         x_grid, y_grid) = resample_on_structutred_grid(data_dict,
                                                        time_point,
                                                        x_min, x_max,
                                                        y_min, y_max,
                                                        **kwargs)
        vector_interpolated.append(quanitity_interpolated)
    return vector_interpolated, x_grid, y_grid


def prepare_mesh(x_grid, y_grid, z_position):
    r"""
    Reshape 2D mesh for 3D vtk writer.
    """
    shape = (x_grid.shape[0], x_grid.shape[1])
    z_grid = np.ones(shape)*z_position
    x_grid = np.expand_dims(np.expand_dims(x_grid, 2), 0)
    y_grid = np.expand_dims(np.expand_dims(y_grid, 2), 0)
    z_grid = np.expand_dims(np.expand_dims(z_grid, 2), 0)
    mesh = np.vstack((x_grid, y_grid, z_grid))
    return mesh


def reshape_vector(vector_x, vector_y, vector_z):
    r"""
    Reshape vector on 2d mesh for 3d vtk writer.
    """
    msg = 'Vector should be on a 2D plane.'
    assert len(vector_x.shape) == len(vector_y.shape) ==  len(vector_z.shape) == 2, msg
    msg = 'Vector components do not have same dimensions.'
    assert vector_x.shape[0] == vector_y.shape[0] == vector_z.shape[0], msg
    assert vector_x.shape[1] == vector_y.shape[1] == vector_z.shape[1], msg
    vector_x = np.expand_dims(np.expand_dims(vector_x, 2), 0)
    vector_y = np.expand_dims(np.expand_dims(vector_y, 2), 0)
    vecotr_z = np.expand_dims(np.expand_dims(vector_z, 2), 0)
    vector = make_obj_array([vector_x, vector_y, vector_z])
    return vector


def reshape_scalar(scalar):
    r"""
    Reshape scalar on 2d mesh for 3d vtk writer.
    """
    msg = 'Scalar should be on a 2D plane.'
    assert len(scalar.shape) == 2, msg
    scalar = np.expand_dims(np.expand_dims(scalar, 2), 0)
    return scalar


def write_to_structured_grid(file_name, data, label, mesh):
    r"""
    Write scalar or vector to a structured grid vtk file.
    """
    write_structured_grid(file_name,
                          mesh,
                          point_data=[(label, data)])