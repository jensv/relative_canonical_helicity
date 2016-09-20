#! /Users/vonderlinden2/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 18:07:14 2016

@author: Jens von der Linden
"""
import argparse
import numpy as np
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator
from scipy.constants import mu_0
from scipy.constants import elementary_charge as q_e
from scipy.constants import proton_mass as m_i
from astropy.convolution import convolve, convolve_fft
from scipy.signal import fftconvolve

from datetime import date
from datetime import datetime
import visit_writer

import structured_3d_vtk as struc_3d
reload(struc_3d)

import os

import ion_current_to_mach_number as ic_to_mach
reload(ic_to_mach)
import read_from_sql


def main(args):
    r"""

    """
    now = datetime.now().strftime("%Y-%m-%d-%H-%M")
    out_dir = '../output/' + now
    try:
        os.makedirs(out_dir)
    except:
        pass

    bx_all_planes = save_idl_quantity_to_unstructured_grids('bx', 'B_x', now,
                                                            *args.bx_extent)
    by_all_planes = save_idl_quantity_to_unstructured_grids('by', 'B_y', now,
                                                            *args.by_extent)
    bz_all_planes = save_idl_quantity_to_unstructured_grids('bz', 'B_z', now,
                                                            *args.bz_extent)
    te_all_planes = save_idl_quantity_to_unstructured_grids('te', 'T_e', now,
                                                            *args.te_extent,
                                                            bounds=args.te_bounds)
    n_all_planes = save_idl_quantity_to_unstructured_grids('n', 'n', now,
                                                           *args.n_extent,
                                                           bounds=args.n_bounds)

    mach_y_all_planes, mach_z_all_planes = prepare_mach_probe_data(args) 
    bx_triangulation, bx_interpolators = give_delaunay_and_interpolator(bx_all_planes)
    by_triangulation, by_interpolators = give_delaunay_and_interpolator(by_all_planes)
    bz_triangulation, bz_interpolators = give_delaunay_and_interpolator(bz_all_planes)
    te_triangulation, te_interpolators = give_delaunay_and_interpolator(te_all_planes)
    n_triangulation, n_interpolators = give_delaunay_and_interpolator(n_all_planes)
    mach_y_triangulation, mach_y_interpolators = give_delaunay_and_interpolator(mach_y_all_planes)
    mach_z_triangulation, mach_z_interpolators = give_delaunay_and_interpolator(mach_z_all_planes)

    (x_min, x_max,
     y_min, y_max,
     z_min, z_max) = args.joint_extent

    mesh = np.meshgrid(np.linspace(x_min, x_max, np.ceil((x_max-x_min)/args.spatial_increment)),
                       np.linspace(y_min, y_max, np.ceil((y_max-y_min)/args.spatial_increment)),
                       np.linspace(z_min, z_max, np.ceil((z_max-z_min)/args.spatial_increment)))

    ## Specific example of velocity calculation
    ##mach_y_interpolator = mach_y_interpolators[70]
    ##mach_z_interpolator = mach_z_interpolators[70]
    ##te_interpolator = te_interpolators[0]
    ##mach_y = scalar_on_mesh(mach_y_interpolator, mesh[:2])
    ##mach_z = scalar_on_mesh(mach_z_interpolator, mesh)
    ##te = scalar_on_mesh(te_interpolator, mesh)
    ##u_i_y = np.sqrt(te*q_e/m_i)*mach_y
    ##u_i_z = np.sqrt(te*q_e/m_i)*mach_z
    ##u_i_y = np.reshape(u_i_y, mesh[0].shape)
    ##u_i_z = np.reshape(u_i_z, mesh[0].shape)
    ##u_i_y = u_i_y[3:, :]
    ##u_i_z = u_i_z[3:, :]


    mesh = np.meshgrid(np.linspace(x_min, x_max, np.ceil((x_max-x_min)/args.spatial_increment)),
                       np.linspace(y_min, y_max, np.ceil((y_max-y_min)/args.spatial_increment)),
                       np.linspace(z_min, z_max, np.ceil((z_max-z_min)/args.spatial_increment)))

    mesh_wo_edges = remove_edges_mesh([np.array(mesh[0]),
                                       np.array(mesh[1]),
                                       np.array(mesh[2])])
    ones = np.ones(mesh_wo_edges[0].shape)

    quantity_names = ['B_x', 'B_y', 'B_z',
                      'B_norm_x', 'B_norm_y', 'B_norm_z',
                      'j_x', 'j_y', 'j_z', 'n', 'Te',
                      'u_i_term1_x', 'u_i_term1_y', 'u_i_term1_z',
                      'u_e_norm_x', 'u_e_norm_y', 'u_e_norm_z',
                      'w_i_term1_x', 'w_i_term1_y', 'w_i_term1_z',
                      'w_i_term2_x', 'w_i_term2_y', 'w_i_term2_z',
                      'u_i_term1_x_constant_density', 'u_i_term1_y_constant_density', 'u_i_term1_z_constant_density',
                      'w_i_term1_x_constant_density', 'w_i_term1_y_constant_density', 'w_i_term1_z_constant_density',
                      'ones']

    for time_point in xrange(len(bx_interpolators)):
        print time_point
        bx_interpolator = bx_interpolators[time_point]
        by_interpolator = by_interpolators[time_point]
        bz_interpolator = bz_interpolators[time_point]
        te_interpolator = te_interpolators[time_point]
        n_interpolator = n_interpolators[time_point]

        bx_derivative = triangulate_derivatives(mesh, bx_triangulation, bx_interpolator,
                                                increment=args.derivative_increment)
        bx_derivative = remove_edges_derivative_meshes(bx_derivative)
        by_derivative = triangulate_derivatives(mesh, by_triangulation, by_interpolator,
                                                increment=args.derivative_increment)
        by_derivative = remove_edges_derivative_meshes(by_derivative)
        bz_derivative = triangulate_derivatives(mesh, bz_triangulation, bz_interpolator,
                                                increment=args.derivative_increment)
        bz_derivative = remove_edges_derivative_meshes(bz_derivative)


        current = current_on_mesh([bx_derivative,
                                   by_derivative,
                                   bz_derivative])
        b_field, b_field_norm = b_field_on_mesh([bx_interpolator,
                                                 by_interpolator,
                                                 bz_interpolator], mesh_wo_edges, bias=args.bias_field_magnitude)

        temperature = scalar_on_mesh(te_interpolator, mesh_wo_edges)
        density = scalar_on_mesh(n_interpolator, mesh_wo_edges)

        current = np.asarray(current)
        density = np.asarray(density)
        b_field_norm = np.asarray(b_field_norm)

        density = boxcar_filter_quantity_mesh(density, args.filter_width)

        for direction in xrange(len(current)):
            current[direction] = boxcar_filter_quantity_mesh(current[direction], args.filter_width)

        density_constant = args.density_constant_factor*np.ones(density.shape)

        ion_velocity_term_1 = calc_ion_velocity_term_1(current, density, q_e)
        ion_velocity_term_1_constant_density = calc_ion_velocity_term_1(current, density_constant, q_e)
        ion_velocity_term_2 = calc_ion_velocity_term_2(b_field_norm, args.alpha)

        ion_vorticity_term_1 = calc_ion_vorticity_term_1(current, density, q_e, mesh_wo_edges)
        ion_vorticity_term_1_constant_density = calc_ion_vorticity_term_1(current, density_constant, q_e, mesh_wo_edges)
        ion_vorticity_term_2 = calc_ion_vorticity_term_2(b_field_norm, args.alpha, mesh_wo_edges)

        for direction in xrange(len(ion_vorticity_term_1)):
            ion_vorticity_term_1[direction] = boxcar_filter_quantity_mesh(ion_vorticity_term_1[direction], args.filter_width)
            ion_vorticity_term_1_constant_density[direction] = boxcar_filter_quantity_mesh(ion_vorticity_term_1_constant_density[direction],
                                                                                           args.filter_width)
            ion_vorticity_term_2[direction] = boxcar_filter_quantity_mesh(ion_vorticity_term_2[direction], args.filter_width)


        fields = (list(b_field) + list(b_field_norm) + list(current) +
                  [density] + [temperature] +
                  list(ion_velocity_term_1) + list(ion_velocity_term_2) +
                  list(ion_vorticity_term_1) + list(ion_vorticity_term_2) +
                  list(ion_velocity_term_1_constant_density) +
                  list(ion_vorticity_term_1_constant_density) +
                  [ones])

        numpy_archive_name = out_dir + args.output_prefix + str(time_point).zfill(4) + '.npz'
        save_to_numpy_mesh(mesh_wo_edges, fields[5:9], quantity_names[5:9], numpy_archive_name)

        x, y, z, variables = prepare_for_rectilinear_grid(mesh_wo_edges, fields,
                                                          quantity_names)

        write_fields_and_currents_to_structured_mesh(now, args.output_prefix, x, y, z, variables, time_point)


def parse_args():
    r"""
    """
    parser = argparse.ArgumentParser(description='Create VTK files of canonical quantities')
    parser.add_argument('--spatial_increment', help='Spatial increment of output file grids', type=float, default=0.001)
    parser.add_argument('--filter_width', help='width of boxcar filter for derivatives', type=float, default=15)
    parser.add_argument('--alpha', help='alpha parameter in ion velocity, usually 1 can be set later by VisIt',
                        type=float, default=1)
    parser.add_argument('--derivative_increment',
                        help='spatial increment used to determine tetrahedron derivative of Delaunay',
                        type=float, default=0.0001)
    parser.add_argument('--density_constant_factor', help='value used for constant density parameters ', type=float,
                        default=1e18)
    parser.add_argument('--bx_extent', help='spatial extent of Bx measurements', nargs=6, type=float,
                        default=[-0.032, 0.028, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--by_extent', help='spatial extent of By measurements', nargs=6, type=float,
                        default=[-0.032, 0.028, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--bz_extent', help='spatial extent of Bz measurements', nargs=6, type=float,
                        default=[-0.032, 0.028, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--te_extent', help='spatial extent of temperature measurements', nargs=6, type=float,
                        default=[-0.026, 0.028, -0.03, 0.028, 0.249, 0.416])
    parser.add_argument('--te_bounds', help='sensible bounds for temperature measurements', nargs=2, type=float,
                        default=[1e-3, 1e3])
    parser.add_argument('--n_extent', help='spatial extent of density measurements', nargs=6, type=float,
                        default=[-0.026, 0.028, -0.03, 0.028, 0.249, 0.416])
    parser.add_argument('--n_bounds', help='sensible bounds for density measurements', nargs=2, type=float,
                        default=[1e3, 1e22])
    parser.add_argument('--joint_extent', help='overlapping spatial extent of all parameters', nargs=6, type=float,
                        default=[-0.022, 0.024, -0.02, 0.018, 0.249, 0.416])
    parser.add_argument('--mach_time_steps', help='# of time steps to extract from one gyration', type=int,
                        default=250)
    parser.add_argument('--shot_database', help='path to shot database',
                        default='/Users/vonderlinden2/rsx_analysis/shots_database/source/shots.db')
    parser.add_argument('--table_name', help='name of sql table',
                        default='Shots')
    parser.add_argument('--min_spectral',
                        help='minimum spectral energy around gyration frequency to include shot', type=float,
                        default=1.6e-8)
    parser.add_argument('--mach_y_extent', help='spatial extent of mach measurements to include', nargs=6, type=float,
                        default=[-0.052, 0.052, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--mach_z_extent', help='spatial extent of mach measurements to include', nargs=6, type=float,
                        default=[-0.032, 0.032, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--mach_bounds', help='bounds on mach measurements', nargs=2, type=float,
                        default=[-10, 10])
    parser.add_argument('--output_prefix', help='prefix of output files', default='Bdot_triple_probe_quantities')
    parser.add_argument('--bias_field_magnitude', help='magnitude of axial bias magnetic field', default=0.02)
    args = parser.parse_args()
    return args


def prepare_mach_probe_data(args):
    r"""
    """
    timesteps = args.mach_time_steps
    database = args.shot_database
    table = args.table_name
    min_spectral_density = args.min_spectral

    z_direction_1, z_direction_2 = 0, 180
    y_direction_1, y_direction_2 = 90, 270
    angle_signs = {0: 1,
                   180: -1,
                   90: -1,
                   0: 1}

    condition_z_0416 = ("campaigns = 'mach_probe_plane_campaign_1'"
                        " AND fiducial_pre_crowbar_gyration_spectral_density > "
                        + str(min_spectral_density) +
                        " AND mach_signals_exist = 1"
                        " AND (mach_orientation = " + str(z_direction_1) +
                        " OR mach_orientation = " + str(z_direction_2) + ")")

    condition_y_0416 = ("campaigns = 'mach_probe_plane_campaign_1'"
                        " AND fiducial_pre_crowbar_gyration_spectral_density > "
                        + str(min_spectral_density) +
                        " AND mach_signals_exist = 1"
                        " AND (mach_orientation = " + str(y_direction_1) +
                        " OR mach_orientation = " + str(y_direction_2) + ")")

    cursor, connection = read_from_sql.cursor_with_rows(condition_z_0416,
                                                        database,
                                                        table)
    z_0416_shots = cursor.fetchall()
    cursor.close()
    connection.close()

    cursor, connection = read_from_sql.cursor_with_rows(condition_y_0416,
                                                        database,
                                                        table)
    y_0416_shots = cursor.fetchall()
    cursor.close()
    connection.close()

    condition_z_302 = ("campaigns = 'mach_probe_plane_campaign_2'"
                       " AND fiducial_pre_crowbar_gyration_spectral_density > "
                       + str(min_spectral_density) +
                       " AND mach_signals_exist = 1"
                       " AND (mach_orientation = " + str(z_direction_1) +
                       " OR mach_orientation = " + str(z_direction_2) + ")")

    cursor, connection = read_from_sql.cursor_with_rows(condition_z_302,
                                                        database,
                                                        table)
    z_0302_shots = cursor.fetchall()
    cursor.close()
    connection.close()

    mach_z_0416_measurements = ic_to_mach.run_mach_analysis(z_0416_shots,
                                                            timesteps,
                                                            angle_signs)
    mach_y_0416_measurements = ic_to_mach.run_mach_analysis(y_0416_shots,
                                                            timesteps,
                                                            angle_signs)
    mach_z_0302_measurements = ic_to_mach.run_mach_analysis(z_0302_shots,
                                                            timesteps,
                                                            angle_signs)

    mach_z_0416_measurements['delays'] = np.arange(timesteps)
    mach_y_0416_measurements['delays'] = np.arange(timesteps)
    mach_z_0302_measurements['delays'] = np.arange(timesteps)

    mach_z_0416_measurements = struc_3d.average_duplicate_points(mach_z_0416_measurements)
    mach_y_0416_measurements = struc_3d.average_duplicate_points(mach_y_0416_measurements)
    mach_z_0302_measurements = struc_3d.average_duplicate_points(mach_z_0302_measurements)

    mach_y_measurements = {0.416: mach_y_0416_measurements}
    mach_z_measurements = {0.302: mach_z_0302_measurements,
                           0.416: mach_z_0416_measurements}

    mach_y_all_planes = save_quantity_to_unstructured_grids(mach_y_measurements,
                                                            'Mach_y', 'Mach_y', '2016-07-26',
                                                            planes=[0.416],
                                                            *args.mach_y_extent,
                                                            bounds=args.mach_bounds)

    mach_z_all_planes = save_quantity_to_unstructured_grids(mach_z_measurements,
                                                            'Mach_z', 'Mach_z', '2016-07-26',
                                                            planes=[0.302, 0.416],
                                                            *args.mach_z_extent,
                                                            bounds=args.mach_bounds)
    mach_y_all_planes = remove_nan_points(mach_y_all_planes)
    mach_z_all_planes = remove_nan_points(mach_z_all_planes)

    return mach_y_all_planes, mach_z_all_planes


def prepare_idl_quantity(name, planes, bounds=None):
    r"""
    """
    measurements = struc_3d.read_idl(name, )
    if bounds:
        measurements = struc_3d.remove_points_out_of_bounds(measurements, 
                                                            bounds[0], bounds[1], planes)
    for plane in planes:
        measurements[plane] = struc_3d.average_duplicate_points(measurements[plane])
    return measurements


def combine_all_planes(measurements, planes):
    r"""
    """
    all_planes = dict(measurements[planes[0]])
    all_planes['z_out'] = planes[0]*np.ones(all_planes['x_out'].size) 
    for plane in planes[1:]:
        for key in ['x_out', 'y_out']:
            all_planes[key] = np.concatenate((all_planes[key], 
                                              measurements[plane][key]))
        all_planes['z_out'] = np.concatenate((all_planes['z_out'],
                                              plane*np.ones(measurements[plane]['x_out'].size)))
        for key in ['std', 'a_out']:
            for time_point in xrange(measurements[plane]['delays'].size):
                all_planes[key][time_point] = np.concatenate((all_planes[key][time_point], 
                                                              measurements[plane][key][time_point]))
    return all_planes


def remove_points_outside_convex(measurements,
                                 x_min= -0.025, x_max=0.025,
                                 y_min=-0.016, y_max=0.017,
                                 z_min=0.249, z_max=0.416):
    r"""
    """
    points = np.dstack((measurements['x_out'], 
                        measurements['y_out'], 
                        measurements['z_out']))[0]
    outside_convex_volume = np.where(np.logical_or.reduce((points[:, 0] < x_min, 
                                                           points[:, 0] > x_max, 
                                                           points[:, 1] < y_min,
                                                           points[:, 1] > y_max,
                                                           points[:, 2] < z_min,
                                                           points[:, 2] > z_max)))[0]
    for key in ['x_out', 'y_out', 'z_out']:
        measurements[key] = np.delete(measurements[key], 
                                      outside_convex_volume)
    for key in ['std', 'a_out']:
        for time_point in xrange(measurements['delays'].size):
            measurements[key][time_point] = np.delete(measurements[key][time_point], 
                                                      outside_convex_volume)
    return measurements


def remove_nan_points(measurements):
    r"""
    """
    nan_positions = np.where(np.isnan(measurements['a_out']))[1]
    for key in ['x_out', 'y_out', 'z_out']:
        measurements[key] = np.delete(measurements[key], 
                                      nan_positions)
    for key in ['std', 'a_out']:
        measurements[key] = np.delete(measurements[key], 
                                      nan_positions, axis=1)
    return measurements


def prepare_for_unstructured_vtk(measurements, quantity_name):
    r"""
    """
    points = np.dstack((measurements['x_out'], 
                        measurements['y_out'], 
                        measurements['z_out']))[0]
    triangulation = Delaunay(points)
    points = tuple(points.ravel())
    connectivity = tuple([(visit_writer.tetrahedron, int(simplex[0]), 
                           int(simplex[1]), int(simplex[2]), int(simplex[3])) 
                           for simplex in triangulation.simplices])
    variables_all_time = []
    for time_point in xrange(measurements['delays'].size):
        variables = ((quantity_name, 1, 1, tuple(measurements['a_out'][time_point])),)
        variables_all_time.append(variables)
    return points, connectivity, variables_all_time


def give_delaunay_and_interpolator(measurements):
    r"""
    """
    if np.unique(measurements['z_out']).size < 2:
        points = np.dstack((measurements['x_out'], 
                            measurements['y_out']))[0]
    else:
        points = np.dstack((measurements['x_out'], 
                            measurements['y_out'], 
                            measurements['z_out']))[0]
    triangulation = Delaunay(points)
    interpolators = []
    for time_point in xrange(measurements['delays'].size):
        interpolators.append(LinearNDInterpolator(points, measurements['a_out'][time_point]))
    return triangulation, interpolators


def write_all_time_unstructured(file_prefix, points, connectivity, variables_all_time):
    r"""
    """
    for time_point in xrange(len(variables_all_time)):
        path = file_prefix + str(time_point).zfill(4)
        visit_writer.WriteUnstructuredMesh(path, 1, points,
                                           connectivity, 
                                           variables[time_point])


def save_idl_quantity_to_unstructured_grids(idl_quantity_name,
                                            visit_quantity_name,
                                            date,
                                            x_min= -0.021, x_max=0.014,
                                            y_min=-0.016, y_max=0.017,
                                            z_min=0.249, z_max=0.416,
                                            planes=[0.249, 0.302, 0.357, 0.416],
                                            file_name_descriptor='_all_planes_convex_unstructured_grid_',
                                            bounds=None):
    r"""
    """
    file_prefix = '../output/' + date + '/' + visit_quantity_name + file_name_descriptor
    measurements = prepare_idl_quantity(idl_quantity_name, planes, bounds=bounds)
    all_planes = combine_all_planes(measurements, planes) 
    all_planes = remove_points_outside_convex(all_planes, 
                                              x_min=x_min, x_max=x_max,
                                              y_min=y_min, y_max=y_max,
                                              z_min=z_min, z_max=z_max)
    assert len(all_planes['x_out']) == len(all_planes['y_out']) == len(all_planes['z_out'])
    #(points, 
    # connectivity, 
    # variables_all_time) = prepare_for_unstructured_vtk(all_planes, visit_quantity_name)
    #assert len(points) == len(all_planes['x_out'])*3
    #assert len(variables_all_time[0][0][3]) == len(all_planes['x_out'])
    #write_all_time_unstructured(file_prefix, points, 
    #                            connectivity, variables_all_time)
    return all_planes 


def save_quantity_to_unstructured_grids(measurements,
                                        idl_quantity_name,
                                        visit_quantity_name,
                                        date,
                                        x_min= -0.021, x_max=0.014,
                                        y_min=-0.016, y_max=0.017,
                                        z_min=0.249, z_max=0.416,
                                        planes=[0.249, 0.302, 0.357, 0.416],
                                        file_name_descriptor='_all_planes_convex_unstructured_grid_',
                                        bounds=None):
    r"""
    """
    file_prefix = '../output/' + date + '/' + visit_quantity_name + file_name_descriptor
    all_planes = combine_all_planes(measurements, planes) 
    all_planes = remove_points_outside_convex(all_planes, 
                                              x_min=x_min, x_max=x_max,
                                              y_min=y_min, y_max=y_max,
                                              z_min=z_min, z_max=z_max)
    assert len(all_planes['x_out']) == len(all_planes['y_out']) == len(all_planes['z_out'])
    #(points, 
    # connectivity, 
    # variables_all_time) = prepare_for_unstructured_vtk(all_planes, visit_quantity_name)
    #assert len(points) == len(all_planes['x_out'])*3
    #assert len(variables_all_time[0][0][3]) == len(all_planes['x_out'])
    #write_all_time_unstructured(file_prefix, points, 
    #                            connectivity, variables_all_time)
    return all_planes 


def triangulate_derivatives(mesh, triangulation, interpolator, increment=0.00001):
    r"""
    """
    d_dx_mesh = np.zeros(mesh[0].shape)
    d_dy_mesh = np.zeros(mesh[0].shape)
    d_dz_mesh = np.zeros(mesh[0].shape)
    points = np.dstack((mesh[0].ravel(),
                        mesh[1].ravel(),
                        mesh[2].ravel()))
    values = interpolator(points)
    d_dx = (interpolator(points + np.asarray([increment, 0, 0])) - values) / increment
    d_dy = (interpolator(points + np.asarray([0, increment, 0])) - values) / increment
    d_dz = (interpolator(points + np.asarray([0, 0, increment])) - values) / increment
    d_dx_mesh = d_dx.reshape(mesh[0].shape)
    d_dy_mesh = d_dy.reshape(mesh[0].shape)
    d_dz_mesh = d_dz.reshape(mesh[0].shape)
    derivative_meshes = ([d_dx_mesh, d_dy_mesh, d_dz_mesh])
    return derivative_meshes


def remove_edges_derivative_meshes(derivative_meshes,
                                   x_start=2, x_end=None, 
                                   y_start=0, y_end=-2, 
                                   z_start=0, z_end=-1):
    r"""
    """
    for index in xrange(len(derivative_meshes)):
        derivative_meshes[index] = derivative_meshes[index][y_start:y_end, x_start:x_end, z_start:z_end]
    return derivative_meshes


def remove_edges_vector_quantity_meshes(quantity_meshes,
                                        x_start=2, x_end=None, 
                                        y_start=0, y_end=-2, 
                                        z_start=0, z_end=-1):
    r"""
    """
    for index in xrange(len(quantity_meshes)):
        quantity_meshes[index] = quantity_meshes[index][y_start:y_end, x_start:x_end, z_start:z_end]
    return quantity_meshes


def remove_edges_scalar_quantity_meshes(quantity_mesh,
                                        x_start=2, x_end=None, 
                                        y_start=0, y_end=-2, 
                                        z_start=0, z_end=-1):
    r"""
    """
    quantity_mesh = quantity_mesh[y_start:y_end, x_start:x_end, z_start:z_end]
    return quantity_mesh


def remove_edges_mesh(mesh,
                      x_start=2, x_end=None, 
                      y_start=0, y_end=-2, 
                      z_start=0, z_end=-1):
    r"""
    """
    for index in xrange(len(mesh)):
        mesh[index] = mesh[index][y_start:y_end, x_start:x_end, z_start:z_end]
    return mesh


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


def save_to_numpy_mesh(mesh, quantities, quantity_names, filename):
    r"""
    """
    dict_to_save = {'mesh_x': mesh[0],
                    'mesh_y': mesh[1],
                    'mesh_z': mesh[2]}
    for i, quanitity in enumerate(quantities):
        dict_to_save[quantity_names[i]] = quantities
    np.savez(filename, **dict_to_save)


def current_on_mesh(derivative_meshes_all_directions, mu_0=mu_0):
    r"""
    """
    dBx_dy = derivative_meshes_all_directions[0][1]
    dBx_dz = derivative_meshes_all_directions[0][2]
    dBy_dx = derivative_meshes_all_directions[1][0]
    dBy_dz = derivative_meshes_all_directions[1][2]
    dBz_dx = derivative_meshes_all_directions[2][0]
    dBz_dy = derivative_meshes_all_directions[2][1]
    j_x = 1./(mu_0) * (dBz_dy - dBy_dz)
    j_y = 1./(mu_0) * (dBx_dz - dBz_dx)
    j_z = 1./(mu_0) * (dBy_dx - dBx_dy)
    current = [j_x, j_y, j_z]    
    return current


def b_field_on_mesh(interpolator_all_directions, mesh, bias=2e-2):
    r"""
    """
    shape = mesh[0].shape
    points = np.swapaxes(np.asarray([mesh[0].ravel(), mesh[1].ravel(), mesh[2].ravel()]), 0, 1)
    b_x = interpolator_all_directions[0](points)  
    b_y = interpolator_all_directions[1](points)
    b_z = interpolator_all_directions[2](points)
    b_x = np.resize(b_x, shape)
    b_y = np.resize(b_y, shape)
    b_z = np.resize(b_z, shape)
    b_z += bias 
    b_x_norm = b_x / np.sqrt(b_x**2. + b_y**2. + b_z**2.)
    b_y_norm = b_y / np.sqrt(b_x**2. + b_y**2. + b_z**2.)
    b_z_norm = b_z / np.sqrt(b_x**2. + b_y**2. + b_z**2.)    
    b_field = [b_x, b_y, b_z]
    b_field_norm = [b_x_norm, b_y_norm, b_z_norm]
    return b_field, b_field_norm


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


def write_fields_and_currents_to_structured_mesh(date, visit_file_name, 
                                                 x, y, z, data, time_point):
    r"""
    """
    file_prefix = '../output/' + date + '/' + visit_file_name
    path = file_prefix + str(time_point).zfill(4)  
    print path
    visit_writer.WriteRectilinearMesh(path, 1, x, y, z, data)


def bdot_probe_extent():
    r"""
    """
    x_min = -0.026
    x_max = 0.024
    y_min = -0.02 
    y_max = 0.028
    z_min = 0.249 
    z_max = 0.416
    return x_min, x_max, y_min, y_max, z_min, z_max

def joint_bdot_tp_extent():
    r"""
    """
    x_min = -0.022
    x_max = 0.024
    y_min = -0.02 
    y_max = 0.024
    z_min = 0.249 
    z_max = 0.416
    return x_min, x_max, y_min, y_max, z_min, z_max

def joint_mach_bdot_tp_extent():
    r"""
    """
    x_min = -0.022
    x_max = 0.024
    y_min = -0.02 
    y_max = 0.018
    z_min = 0.249 
    z_max = 0.416
    return x_min, x_max, y_min, y_max, z_min, z_max


def calc_ion_velocity_term_1(current, density, charge):
    r"""
    """
    denominator = density*charge
    term = [current[0]/denominator,
            current[1]/denominator,
            current[2]/denominator]
    return term


def calc_ion_velocity_term_2(b_field_norm, alpha):
    r"""
    """
    return alpha*b_field_norm


def curl_on_mesh(quantity, mesh):
    r"""
    """
    dx = mesh[0][0, 1, 0] - mesh[0][0, 0, 0]
    dy = mesh[1][1, 0, 0] - mesh[1][0, 0, 0]
    dz = mesh[2][0, 0, 1] - mesh[2][0, 0, 0]
    dx_dy = np.gradient(quantity[0], axis=1)/dy
    dx_dz = np.gradient(quantity[0], axis=2)/dz
    dy_dx = np.gradient(quantity[1], axis=0)/dx
    dy_dz = np.gradient(quantity[1], axis=2)/dz
    dz_dx = np.gradient(quantity[2], axis=0)/dx
    dz_dy = np.gradient(quantity[2], axis=1)/dy
    curl_x = dz_dy - dy_dz
    curl_y = dx_dz - dz_dx
    curl_z = dy_dx - dx_dy
    return [curl_x, curl_y, curl_z]


def calc_ion_vorticity_term_1(current, density, charge, mesh, filt=None):
    r"""
    """
    return curl_on_mesh(calc_ion_velocity_term_1(current, density, charge), mesh)


def calc_ion_vorticity_term_2(b_field_norm, alpha, mesh):
    r"""
    """
    return curl_on_mesh(calc_ion_velocity_term_2(b_field_norm, alpha), mesh)

def boxcar_filter_quantity_mesh(quantity, width):
    r"""
    """
    nan_indexes = np.where(np.isnan(quantity))
    quantity[nan_indexes] = 0
    boxcar = np.ones((width, width, width))
    boxcar /= boxcar.sum()
    return fftconvolve(quantity, boxcar, mode='same')


def boxcar_filter_quantity_mesh(quantity, width):
    r"""
    """
    nan_indexes = np.where(np.isnan(quantity))
    quantity[nan_indexes] = 0
    boxcar = np.ones((width, width, width))
    boxcar /= boxcar.sum()
    return fftconvolve(quantity, boxcar, mode='same')


if __name__ == '__main__':
    args = parse_args()
    main(args) 
