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
import sys


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
                                 x_points=25,
                                 y_points=25,
                                 method='cubic'):
    r"""
    Resample quantity from measurement grid to structured grid.
    """
    x_points = np.linspace(x_min, x_max, x_points)
    y_points = np.linspace(y_min, y_max, y_points)
    x_grid, y_grid = np.meshgrid(x_points, y_points)
    quantity_interpolated = griddata(np.dstack((data_dict['x_out'],
                                                data_dict['y_out']))[0],
                                               data_dict['a_out'][time_point],
                                               (x_grid, y_grid),
                                               method=method)
    quantity_interpolated = quantity_interpolated[x_slice, y_slice]
    x_grid = x_grid[x_slice, y_slice]
    y_grid = y_grid[x_slice, y_slice]
    return quantity_interpolated, x_grid, y_grid


def remove_nans(quantity, x_grid, y_grid):
    r"""
    """
    nan_positions = np.isnan(quantity)
    if np.sum(nan_positions) == 0:
        return quantity, x_grid, y_grid

    while not np.sum(np.isnan(quantity)) == 0:
        shape = quantity.shape
        nan_to_remove = np.isnan(quantity)[0]
        rows_to_remove = np.max([nan_to_remove[0] + 1, shape[0] - nan_to_remove[0]])
        colums_to_remove = np.max([nan_to_remove[1] + 1, shape[1] - nan_to_remove[1]])
        axis = 0 if rows_to_remove >= columns_to_remove else 1
        if nan_to_remove[axis] + 1 >= shape[axis] - nan_to_remove[axis]:
            indexes = np.s_[nan_to_remove[axis] + 1:]
        else:
            indexes = np.s_[:-(shape[axis] - nan_to_remove[axis])]
        np.delete(quantity, indexes, axis=axis)
        np.delete(x_grid, indexes, axis=axis)
        np.delete(y_grid, indexes, axis=axis)
    return quantity, x_grid, y_grid


def remove_nans_vector(vector, x_grid, y_grid):
    r"""
    """
    nan_positions = np.isnan(vector)
    if np.sum(nan_positions) == 0:
        return vector, x_grid, y_grid
    directions = set([0, 1, 2])
    for direction in range(3):
        other_direc = list(directions.difference([direction]))
        quantity = vector[direction]
        other_quantities = [vector[other_direc[0]], vector[other_direc[1]]]
        while not np.sum(np.isnan(quantity)) == 0:
            shape = quantity.shape
            nan_to_remove = np.isnan(quantity)[0]
            rows_to_remove = np.max([nan_to_remove[0] + 1, shape[0] - nan_to_remove[0]])
            colums_to_remove = np.max([nan_to_remove[1] + 1, shape[1] - nan_to_remove[1]])
            axis = 0 if rows_to_remove >= columns_to_remove else 1
            if nan_to_remove[axis] + 1 >= shape[axis] - nan_to_remove[axis]:
                indexes = np.s_[nan_to_remove[axis] + 1:]
            else:
                indexes = np.s_[:-(shape[axis] - nan_to_remove[axis])]
            np.delete(quantity, indexes, axis=axis)
            np.delete(other_quantities[0], indexes, axis=axis)
            np.delete(other_quantities[1], indexes, axis=axis)
            np.delete(x_grid, indexes, axis=axis)
            np.delete(y_grid, indexes, axis=axis)
        vector[direction] = quantity
        vector[other_direc[0]] = other_quantities[0]
        vector[other_direc[1]] = other_quantities[1]
    return vector, x_grid, y_grid


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


def main():
    r"""
    Called from command line. Create a vector or scalar vtk file from idl sav
    files.
    """
    input_dict = handle_args()
    if input_dict['data_type'] == 'vector':
        measurements_x = read_data(input_dict['x_input_path'])
        measurements_y = read_data(input_dict['y_input_path'])
        measurements_z = read_data(input_dict['z_input_path'])
        vector_dicts = [measurements_x, measurements_y, measurements_z]
        (x_min, x_max, y_min, y_max) = determine_sample_bounds(vector_dicts)
        for time_point in input_dict['time_points']:
            quantity_resampled, x_grid, y_grid = resample_vector(vector_dicts[0],
                                                                 vector_dicts[1],
                                                                 vector_dicts[2],
                                                                 time_point,
                                                                 x_min, x_max, y_min, y_max,
                                                                 to_clip=input_dict['to_clip'])
            mesh = prepare_mesh(x_grid, y_grid, input_dict['z_position'])
            vector = reshape_vector(quantity_resampled)
            output_path = input_dict['output_path'] + str(time_point).zfill(6)
            write_to_structured_grid(output_path, vector,
                                     input_dict['symbol'], mesh)
    if input_dict['data_type'] == 'scalar':
        measurements = read_data['input_path']
        (x_min, x_max, y_min, y_max) = determine_sample_bounds(measurements)
        for time_point in input_dict['time_points']:
            quantity_resampled, x_grid, y_grid = resample_scalar(measurements,
                                                                 time_point,
                                                                 x_min, x_max, y_min, y_max,
                                                                 to_clip=input_dict['to_clip'])
            mesh = prepare_mesh(x_grid, y_grid, input_dict['z_position'])
            scalar = reshape_scalar(quantity_resampled)
            output_path = input_dict['output_path'] + str(timepoint).zfill(6)
            write_to_structured_grid(output_path, scalar, input_dict['symbol'],
                                     output_path)


def handle_args():
    r"""
    Put command line args in dictionary.
    """
    data_type = sys.argv[1]
    input_dict = {}
    if data_type == 'vector':
        msg = 'vector type requires 6 arguments: 3 input files, z_position, time_points, output_file, symbol, to_clip'
        assert len(sys.argv) == 10, msg
        for i, key in enumerate(('data_type', 'x_input_path', 'y_input_path',
                                 'z_input_path', 'z_position', 'time_points',
                                 'output_path', 'symbol', 'to_clip')):
            input_dict[key] = sys.argv[1:][i]
    if data_type == 'scalar':
        assert len(sys.argv) == 8, msg
        msg = 'scalar type requires 4 arguments: input file, z_position, time_points, output_file, symbol, to_clip'
        for i, key in enumerate(('data_type', 'z_position', 'input_path',
                                 'time_points', 'output_path', 'symbol',
                                 'to_clip')):
            input_dict[key] = sys.argv[1:][i]
    input_dict['time_points'] = np.arange(input_dict['time_points'])
    input_dict['to_clip'] = [0, input_dcit['to_clip'], 0, input_dict['to_clip']]
    return input_dict


if __name__ == '__main__':
    main()
