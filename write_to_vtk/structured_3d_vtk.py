# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 13:05:43 2016

@author: Jens von der Linden
"""

import numpy as np
from scipy.interpolate import LinearNDInterpolator
import visit_writer

def get_interpolator(points, values):
    r"""
    """
    return LinearNDInterpolator(points, values)


def triangulate_grad(mesh, interpolator, increment=0.00001):
    r"""
    """

    d_dx_mesh = np.zeros(mesh[0].shape)
    d_dy_mesh = np.zeros(mesh[0].shape)
    if len(mesh[0].shape) == 3:
        d_dz_mesh = np.zeros(mesh[0].shape)
        points = np.dstack((mesh[0].ravel(),
                            mesh[1].ravel(),
                            mesh[2].ravel()))
    else:
        points = np.dstack((mesh[0].ravel(),
                            mesh[1].ravel()))
    values = interpolator(points)
    if len(mesh[0].shape) == 3:
        d_dx = (interpolator(points + np.asarray([increment, 0, 0])) - values) / increment
        d_dy = (interpolator(points + np.asarray([0, increment, 0])) - values) / increment
        d_dz = (interpolator(points + np.asarray([0, 0, increment])) - values) / increment
    else:
        d_dx = (interpolator(points + np.asarray([increment, 0])) - values) / increment
        d_dy = (interpolator(points + np.asarray([0, increment])) - values) / increment
    d_dx_mesh = d_dx.reshape(mesh[0].shape)
    d_dy_mesh = d_dy.reshape(mesh[0].shape)
    if len(mesh[0].shape) == 3:
        d_dz_mesh = d_dz.reshape(mesh[0].shape)
        derivative_meshes = ([d_dx_mesh, d_dy_mesh, d_dz_mesh])
    else:
        derivative_meshes = ([d_dx_mesh, d_dy_mesh])
    return derivative_meshes


def add_vacuum_field(field, vacuum_field=0.02):
    r"""
    Add constant axial vacuum_field to vector field.
    """
    field[2] = field[2] + vacuum_field
    return field


def scalar_on_mesh(interpolator, mesh):
    r"""
    """
    shape = mesh[0].shape
    assert (len(mesh) == 3 or len(mesh) == 2)
    if len(mesh) == 3:
        points = np.swapaxes(np.asarray([mesh[0].ravel(), mesh[1].ravel(), mesh[2].ravel()]), 0, 1)
    else:
        points = np.swapaxes(np.asarray([mesh[0].ravel(), mesh[1].ravel()]), 0, 1)
    scalar_all_times = []
    scalar = interpolator(points)
    scalar = np.resize(scalar, shape)
    return scalar


def vector_on_mesh(interpolators, mesh):
    r"""
    """
    vec_x = scalar_on_mesh(interpolators[0], mesh)
    vec_y = scalar_on_mesh(interpolators[1], mesh)
    vec_z = scalar_on_mesh(interpolators[2], mesh)
    return (vec_x, vec_y, vec_z)


def prepare_for_rectilinear_grid(mesh, quantities, quantity_names):
    r"""
    """
    x = tuple(np.unique(mesh[0]))
    y = tuple(np.unique(mesh[1]))
    z = tuple(np.unique(mesh[2]))
    variables = [(quantity_names[index], 1, 1,
                  tuple(np.swapaxes(np.swapaxes(quantities[index], 1, 2), 0, 1).ravel()))
                 for index in xrange(len(quantities))]
    return x, y, z, variables


def write_fields_to_rectilinear_grid(visit_file_prefix,
                                     x, y, z, data, time_point):
    r"""
    """
    path = visit_file_prefix + str(time_point).zfill(4)
    visit_writer.WriteRectilinearMesh(path, 1, x, y, z, data)
