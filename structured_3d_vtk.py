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


def read_idl(quantity, data_path='../../comprehensive_3d_plot/output/2016-04-07/'):
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


def interpolate_vector(grid_points, vector_points, vector_values):
    r"""
    Linearly interpolate vector measurements at grid points.
    """
    interpolated_data_x = griddata(vector_points[0], vector_values[0], grid_points)
    interpolated_data_y = griddata(vector_points[1], vector_values[1], grid_points)
    interpolated_data_z = griddata(vector_points[2], vector_values[2], grid_points)
    return [interpolated_data_x, interpolated_data_y, interpolated_data_z]


def interpolate_scalar(grid_points, scalar_points, scalar_values):
    r"""
    Linearly interpolate scalar measurements at grid points.
    """
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
