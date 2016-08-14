# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 13:05:43 2016

@author: Jens von der Linden
"""

import numpy as np
from pyvisfile.vtk import write_structured_grid
from pytools.obj_array import make_obj_array
import scipy.io.idl as idl
from scipy.interpolate import griddata
from collections import MutableSequence

import sys
sys.path.append('/Users/vonderlinden2/rsx_analysis/read_from_sql')
import read_from_sql
sys.path.append('/Users/vonderlinden2/rsx_analysis/mach_probe_analysis')
import ion_current_to_mach_number as ic_to_mach
sys.path.append('/Users/vonderlinden2/rsx_analysis/time_alignment/source')
import absolute_times as at


def read_idl(quantity, data_path='../../comprehensive_3d_plot/output/2016-08-12/'):
    r"""
    Read idl files for all planes.
    """
    idl_ending = '.sav'
    z249_file = data_path + quantity + '_z249' + idl_ending
    z302_file = data_path + quantity + '_z302' + idl_ending
    z357_file = data_path + quantity + '_z357' + idl_ending
    z416_file = data_path + quantity + '_z416' + idl_ending

    z249_measurements = idl.readsav(z249_file)
    z302_measurements = idl.readsav(z302_file)
    z357_measurements = idl.readsav(z357_file)
    z416_measurements = idl.readsav(z416_file)

    measurements = {0.249: z249_measurements,
                    0.302: z302_measurements,
                    0.357: z357_measurements,
                    0.416: z416_measurements}

    return measurements


def read_points_from_measurement_dict(measurement_dict, time_point, z_planes):
    r"""
    Collectand return measurement points and values from dictionary.
    """
    x_points = np.empty((0))
    y_points = np.empty((0))
    z_points = np.empty((0))
    values = np.empty((0))
    for z_plane in z_planes:
        plane_measurements = measurement_dict[z_plane]
        x_points = np.append(x_points, plane_measurements['x_out'])
        y_points = np.append(y_points, plane_measurements['y_out'])
        z_points = np.append(z_points, np.ones(plane_measurements['x_out'].size)*z_plane)
        values = np.append(values, plane_measurements['a_out'][time_point])

    points = [x_points, y_points, z_points]
    points = np.asarray(points)
    points = np.swapaxes(points, 0, 1)
    return points, values


def average_duplicate_points(data_dict):
    r"""
    Find duplicate points, average them and record standard deviation.
    """
    data_dict['x_out'] = data_dict['x_out'].astype('float64')
    data_dict['y_out'] = data_dict['y_out'].astype('float64')
    data_dict['a_out'] = data_dict['a_out'].astype('float64')
    time_points = data_dict['a_out'].shape[0]
    data = {}
    for idx in xrange(data_dict['x_out'].size):
        location = (data_dict['x_out'][idx], data_dict['y_out'][idx])
        if location in data.keys():
            data[location] = np.column_stack((data[location], data_dict['a_out'][:, idx]))
        else:
            data[location] = data_dict['a_out'][:, idx]

    unique_data_dict = {'x_out': [],
                        'y_out': [],
                        'a_out': [],
                        'std': []}
    for location in data.keys():
        if data[location][0].size > 1:
            unique_data_dict['std'].append(data[location].std(axis=1, ddof=1))
            unique_data_dict['a_out'].append(data[location].mean(axis=1))
        else:
            unique_data_dict['std'].append(np.zeros(time_points))
            unique_data_dict['a_out'].append(data[location])
        unique_data_dict['x_out'].append(location[0])
        unique_data_dict['y_out'].append(location[1])

    unique_data_dict['x_out'] = np.asarray(unique_data_dict['x_out'])
    unique_data_dict['y_out'] = np.asarray(unique_data_dict['y_out'])
    unique_data_dict['a_out'] = np.hsplit(np.asarray(unique_data_dict['a_out']), time_points)
    unique_data_dict['std'] = np.hsplit(np.asarray(unique_data_dict['std']), time_points)
    unique_data_dict['delays'] = data_dict['delays']
    return unique_data_dict


def remove_points_out_of_bounds(data_dict, lower, upper, planes):
    r"""
    Find points out of bounds and remove them.
    """
    new_data_dict = {}
    for plane in planes:
        new_data_dict[plane] = {'x_out': [],
                                'y_out': [],
                                'a_out': []}
        to_remove = []
        for time in xrange(len(data_dict[plane]['a_out'])):
            indexes = (np.where(np.logical_or(data_dict[plane]['a_out'][time] < lower,
                                              data_dict[plane]['a_out'][time] > upper))[0])
            if indexes.size > 0:
                to_remove.append(indexes)
        if len(to_remove) > 0:
            to_remove = np.concatenate(to_remove)
        to_remove = np.unique(to_remove)
        new_data_dict[plane]['x_out'] = np.delete(data_dict[plane]['x_out'], to_remove)
        new_data_dict[plane]['y_out'] = np.delete(data_dict[plane]['y_out'], to_remove)
        for time in xrange(len(data_dict[plane]['a_out'])):
            new_data_dict[plane]['a_out'].append(np.delete(data_dict[plane]['a_out'][time],
                                                 to_remove))
        new_data_dict[plane]['a_out'] = np.asarray(new_data_dict[plane]['a_out'])
        new_data_dict[plane]['delays'] = data_dict[plane]['delays']
    return new_data_dict


def bounded_grid(bounds, spatial_increment=0.003):
    r"""
    Return rectinlinear grid bounded by bounds with spatial_increment.
    """
    (x_min, x_max), (y_min, y_max), (z_min, z_max) = bounds

    x_coord = np.linspace(x_min, x_max, np.ceil((x_max-x_min)/spatial_increment))
    if x_coord.size == 0:
        x_coord = np.asarray([x_min])
    y_coord = np.linspace(y_min, y_max, np.ceil((y_max-y_min)/spatial_increment))
    if y_coord.size == 0:
        y_coord = np.asarray([y_min])
    z_coord = np.linspace(z_min, z_max, np.ceil((z_max-z_min)/spatial_increment))
    if z_coord.size == 0:
        z_coord = np.asarray([z_min])
    sizes = map(np.size, [x_coord, y_coord, z_coord])

    mesh = np.meshgrid(x_coord, y_coord, z_coord, indexing='ij')

    grid_points = np.dstack(map(np.ravel, mesh))[0]

    return grid_points, sizes


def interpolate_vector(grid_points, vector_points, vector_values,
                       method='linear'):
    r"""
    Linearly interpolate vector measurements at grid points.
    """
    interpolated_data_x = griddata(vector_points[0], vector_values[0],
                                   grid_points, method=method)
    interpolated_data_y = griddata(vector_points[1], vector_values[1],
                                   grid_points, method=method)
    interpolated_data_z = griddata(vector_points[2], vector_values[2],
                                   grid_points, method=method)
    return [interpolated_data_x, interpolated_data_y, interpolated_data_z]


def interpolate_scalar(grid_points, scalar_points, scalar_values):
    r"""
    Linearly interpolate scalar measurements at grid points.
    """
    #print 'grid' , grid_points
    #print 'm points' , scalar_points
    #print 'values' , scalar_values
    interpolated_data = griddata(scalar_points, scalar_values, grid_points)
    return interpolated_data


def add_vacuum_field(field, vacuum_field=0.02):
    r"""
    Add constant axial vacuum_field to vector field.
    """
    field[2] = field[2] + vacuum_field
    return field


def prepare_mesh(grid_points, sizes):
    r"""
    Return grid reshaped for vtk writer.
    """
    return grid_points.swapaxes(0,1).reshape((3, sizes[0], sizes[1], sizes[2]))


def prepare_vector(vector, sizes):
    r"""
    Return vector rehsaped for vtk writer.
    """
    vtk_vector_x = np.resize(vector[0], sizes)
    vtk_vector_y = np.resize(vector[1], sizes)
    vtk_vector_z = np.resize(vector[2], sizes)

    vtk_vector_x = np.expand_dims(vtk_vector_x, 0)
    vtk_vector_y = np.expand_dims(vtk_vector_y, 0)
    vtk_vector_z = np.expand_dims(vtk_vector_z, 0)

    vtk_vector = make_obj_array([vtk_vector_x, vtk_vector_y, vtk_vector_z])
    return vtk_vector


def prepare_scalar(scalar, sizes):
    r"""
    Return scalar reshaped for vtk writer.
    """
    vtk_scalar = np.resize(scalar, sizes)
    vtk_scalar = np.expand_dims(vtk_scalar, 0)
    return vtk_scalar


def build_vtk_vector(direction_measurements, full_vtk_grid=None, indices=None,
                     bounds=None, spatial_increment=None,
                     z_planes=[0.249, 0.302, 0.357, 0.416], time_points=21):
    r"""
    Return interpolated vector on a vtk grid.
    """
    mode = ''
    msg = 'Either set full_vtk_grid and indcies or bounds and increment not all 4.'
    assert (((not full_vtk_grid == None and not indices == None) or
             (not bounds == None and not spatial_increment == None)) and
            not (not full_vtk_grid == None and not indices == None and
                 not bounds == None and not spatial_increment == None )), msg
    if not full_vtk_grid == None and not indices == None:
        mode = 'sub_grid'
    if not bounds == None and not spatial_increment == None:
        mode = 'full_grid'

    vtk_vectors = []
    for time_point in xrange(time_points):
        print 'time_point %i' % time_point
        points = []
        values = []
        for measurements in direction_measurements:
            (points_direction,
             values_direction) = read_points_from_measurement_dict(measurements,
                                                                   time_point,
                                                                   z_planes)
            points.append(points_direction)
            values.append(values_direction)

        if mode == 'full_grid':
            grid_points, sizes = bounded_grid(bounds, spatial_increment)
        elif mode == 'sub_grid':
            grid_points, sizes = build_sub_grid_points(full_vtk_grid, indices)
        interpolated_vector = interpolate_vector(grid_points, points, values)
        interpolated_vector = add_vacuum_field(interpolated_vector)

        assert np.sum(np.isnan(interpolated_vector[0])) == 0
        assert np.sum(np.isnan(interpolated_vector[1])) == 0
        assert np.sum(np.isnan(interpolated_vector[2])) == 0

        vtk_grid = prepare_mesh(grid_points, sizes)
        vtk_vector = prepare_vector(interpolated_vector, sizes)

        vtk_vectors.append(vtk_vector)
    vtk_vectors = np.asarray(vtk_vectors)
    return vtk_grid, vtk_vectors


def build_vtk_scalar(measurements, full_vtk_grid=None, indices=None,
                     bounds=None, spatial_increment=None,
                     z_planes=[0.249, 0.302, 0.357, 0.416], time_points=21):
    r"""
    Return interpolated scalar on a vtk grid.
    """
    mode = ''
    msg = 'Either set full_vtk_grid and indcies or bounds and increment not all 4.'
    assert (((not full_vtk_grid == None and not indices == None) or
             (not bounds == None and not spatial_increment == None)) and
            not (not full_vtk_grid == None and not indices == None and
                 not bounds == None and not spatial_increment == None )), msg
    if not full_vtk_grid == None and not indices == None:
        mode = 'sub_grid'
    if not bounds == None and not spatial_increment == None:
        mode = 'full_grid'


    vtk_scalars = []
    for time_point in xrange(time_points):
        print 'time_point %i' % time_point
        (scalar_points,
         scalar_values) = read_points_from_measurement_dict(measurements,
                                                            time_point,
                                                            z_planes)
        to_remove = []
        for i, value in enumerate(scalar_values):
            if np.isnan(value):
                to_remove.append(i)
        scalar_values = np.delete(scalar_values, to_remove)
        scalar_points = np.delete(scalar_points, to_remove, 0)

        if mode == 'full_grid':
            grid_points, sizes = bounded_grid(bounds, spatial_increment)
        elif mode == 'sub_grid':
            grid_points, sizes = build_sub_grid_points(full_vtk_grid, indices)

        if len(z_planes) == 1:
                grid_points = np.delete(grid_points, 2, 1)
                scalar_points = np.delete(scalar_points, 2, 1)
        interpolated_scalar = interpolate_scalar(grid_points, scalar_points,
                                                 scalar_values)

        assert np.sum(np.isnan(interpolated_scalar)) == 0

        if len(z_planes) == 1:
            grid_points = np.insert(grid_points, 2,
                                    np.ones((grid_points.shape[0]))*z_planes[0],
                                    axis=1)

        vtk_grid = prepare_mesh(grid_points, sizes)
        vtk_scalar = prepare_scalar(interpolated_scalar, sizes)

        vtk_scalars.append(vtk_scalar)
    vtk_scalars = np.asarray(vtk_scalars)
    return vtk_grid, vtk_scalars


def build_sub_grid_indices(full_vtk_grid, bounds):
    r"""
    """
    indices = np.argwhere(np.logical_and.reduce((full_vtk_grid[0] >= bounds[0][0],
                                                  full_vtk_grid[0] <= bounds[0][1],
                                                  full_vtk_grid[1] >= bounds[1][0],
                                                  full_vtk_grid[1] <= bounds[1][1],
                                                  full_vtk_grid[2] >= bounds[2][0],
                                                  full_vtk_grid[2] <= bounds[2][1])))
    return indices


def build_vtk_sub_indices(indices):
    r"""
    """
    sizes = (np.unique(indices[:, 0]).size,
             np.unique(indices[:, 1]).size,
             np.unique(indices[:, 2]).size)
    indicies = np.swapaxes(indices, 0, 1)
    sub_grid_indices = [np.resize(indicies[i], sizes) for i in xrange(3)]
    return sub_grid_indices


def build_sub_grid_points(full_vtk_grid, indices):
    r"""
    """
    points = [full_vtk_grid[i][indices[:,0],
                               indices[:,1],
                               indices[:,2]] for i in range(3)]
    points = np.asarray(points)
    points = np.swapaxes(points, 0, 1)
    sizes = (np.unique(indices[:, 0]).size,
             np.unique(indices[:, 1]).size,
             np.unique(indices[:, 2]).size)
    return points, sizes