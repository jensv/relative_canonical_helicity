#! /Users/vonderlinden2/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on Apr 13 2017

@author: Jens von der Linden

Filter RSX measurements:
1)interpolate unstructured grid of measurements to rectilinear grid
2)apply Gaussian filter
3)Resample to unstructured grid
"""
import argparse
import numpy as np
from datetime import date
from datetime import datetime
import os

from scipy.interpolate import LinearNDInterpolator
from scipy import ndimage

from write_to_vtk.read_unstructured_vtk import read_unstructured_vtk
from mach_probe_analysis import ion_current_to_mach_number as ic_to_mach
from read_from_sql import read_from_sql
from write_to_vtk import structured_3d_vtk as struc_3d
from write_to_vtk import prepare_measurements as pm
from write_to_vtk import unstructured_grid as ug


def main(args):
    r"""
    """
    now = datetime.now().strftime("%Y-%m-%d-%H-%M")
    out_dir = '../output/filtered_unstructured_measurements/' + now
    try:
        os.makedirs(out_dir)
    except:
        pass

    planes = [0.249, 0.302, 0.357, 0.416]
    bx_measurements = pm.read_idl('bx')
    by_measurements = pm.read_idl('by')
    bz_measurements = pm.read_idl('bz')
    te_measurements = pm.read_idl('te')
    n_measurements = pm.read_idl('n')
    mach_y_measurements, mach_z_measurements = pm.read_mach_probe_data(args)

    if args.bxby_only:
        bx_all_planes = pm.cut_and_average_quantity(bx_measurements,
                                                    args.bxby_extent, planes)
        by_all_planes = pm.cut_and_average_quantity(by_measurements,
                                                    args.bxby_extent, planes)
    else:
        bx_all_planes = pm.cut_and_average_quantity(bx_measurements,
                                                    args.bx_extent, planes)
        by_all_planes = pm.cut_and_average_quantity(by_measurements,
                                                    args.by_extent, planes)
    bz_all_planes = pm.cut_and_average_quantity(bz_measurements,
                                                args.bz_extent, planes)
    n_all_planes = pm.cut_and_average_quantity(n_measurements,
                                               args.n_extent,
                                               planes,
                                               bounds=args.n_bounds)
    te_all_planes = pm.cut_and_average_quantity(te_measurements, args.te_extent,
                                                planes, bounds=args.te_bounds)
    mach_y_plane = pm.cut_and_average_quantity(mach_y_measurements, args.mach_y_extent,
                                               [0.416], bounds=args.mach_bounds)
    mach_z_plane = pm.cut_and_average_quantity(mach_z_measurements, args.mach_z_extent,
                                               [0.416], bounds=args.mach_bounds)


    n_three_planes = pm.remove_plane(0.302, n_all_planes)
    te_three_planes = pm.remove_plane(0.302, te_all_planes)

    if args.bxby_only:
        bxby_grid_x_points = int((args.bxby_filter_extent[1] - args.bxby_filter_extent[0])/
                                  args.filter_spatial_increment)
        bxby_grid_y_points = int((args.bxby_filter_extent[3] - args.bxby_filter_extent[2])/
                                  args.filter_spatial_increment)
        bxby_grid = np.meshgrid(np.linspace(args.bxby_filter_extent[0],
                                            args.bxby_filter_extent[1],
                                            bxby_grid_x_points),
                                np.linspace(args.bxby_filter_extent[2],
                                            args.bxby_filter_extent[3],
                                            bxby_grid_y_points))
        single_plane = np.unique(bx_all_planes['z_out'])[0]
        bx_filtered = filter_unstructured_data(bxby_grid, bx_all_planes,
                                               single_plane=single_plane,
                                               filter_sigma=args.filter_sigma,
                                               filter_truncate=args.filter_truncate)

        by_filtered = filter_unstructured_data(bxby_grid, by_all_planes,
                                               single_plane=single_plane,
                                               filter_sigma=args.filter_sigma,
                                               filter_truncate=args.filter_truncate)

    ug.save_to_unstructured_grid(bx_filtered, 'bx', out_dir,
                                 prefix=args.output_prefix)
    ug.save_to_unstructured_grid(by_filtered, 'by', out_dir,
                                 prefix=args.output_prefix)
    #ug.save_to_unstructured_grid(bz_all_planes, 'bz', out_dir,
    #                             prefix=args.output_prefix)
    #ug.save_to_unstructured_grid(te_three_planes, 'te', out_dir,
    #                             prefix=args.output_prefix)
    #ug.save_to_unstructured_grid(n_three_planes, 'n', out_dir,
    #                             prefix=args.output_prefix)
    #ug.save_to_unstructured_grid(mach_y_plane, 'mach_y', out_dir,
    #                             prefix=args.output_prefix)
    #ug.save_to_unstructured_grid(mach_z_plane, 'mach_z', out_dir,
    #                             prefix=args.output_prefix)


def filter_unstructured_data(grid, measurements, filter_sigma=None,
                             single_plane=None, filter_truncate=None):
    r"""
    Filter data on unstructured grid.

    Interpolate data to rectilinear grid, filter, resample
    onto unstructured grid.
    """
    (planes, points_by_plane,
     values_by_plane) = extract_planes(measurements)
    if single_plane:
        planes = [single_plane]
    filtered_values_by_plane = []
    delays = measurements['delays']
    for i, plane in enumerate(planes):
        filtered_by_time_point = interpolate_and_filter_data(points_by_plane[i],
                                                             values_by_plane[i],
                                                             grid, delays,
                                                             filter_sigma=filter_sigma,
                                                             filter_truncate=filter_truncate)
        (points,
         values_by_time_point) = resample_to_unstructured_grid(grid,
                                                               filtered_by_time_point,
                                                               points_by_plane[i], delays)
        filtered_values_by_plane.append(values_by_time_point)
    filtered_measurements = recombine_planes(planes, points_by_plane,
                                             filtered_values_by_plane,
                                             delays)
    return filtered_measurements


def extract_planes(measurements):
    r"""
    Extract measurement points and values by plane from
    measurement dictionaries.
    """
    planes = np.unique(measurements['z_out'])
    points_by_plane = []
    values_by_plane = []
    measurements['a_out'] = np.asarray(measurements['a_out'])
    for plane in planes:
        indexes = np.where(measurements['z_out'] == plane)[0]
        points = np.stack((measurements['x_out'][indexes],
                           measurements['y_out'][indexes]), axis=1)
        values = measurements['a_out'][:, indexes]
        points_by_plane.append(points)
        values_by_plane.append(values)
    return planes, points_by_plane, values_by_plane


def interpolate_and_filter_data(points, values, grid, delays, filter_sigma=None,
                                filter_truncate=None):
    r"""
    Interpolate and filter (with Gaussian)
    """
    print 'grid in interpolate', grid
    filtered_by_time_point = []
    for time_point in xrange(delays.size):
        print 'filter', time_point
        interpolator = struc_3d.get_interpolator(points, values[time_point])
        data = interpolator(grid[0], grid[1])
        #print data.size, data.shape
        #print np.sum(np.isnan(data))
        #print 'nan x', np.unique(grid[0][np.isnan(data)])
        #print 'nan y', np.unique(grid[1][np.isnan(data)])
        #assert np.sum(np.isnan(data)) == 0, 'interpolated data contains nans'
        if filter_sigma:
            if filter_truncate:
                filtered = ndimage.gaussian_filter(data, filter_sigma,
                                                   truncate=filter_truncate)
            else:
                filtered = ndimage.gaussian_filter(data, filter_sigma)
        else:
            filtered = data
        filtered_by_time_point.append(filtered)
    return filtered_by_time_point


def resample_to_unstructured_grid(grid, data, points, delays):
    r"""
    Resample filtered data back to measurement grid.
    """
    values_by_time_point = []
    grid_points_x = grid[0].ravel()
    grid_points_y = grid[1].ravel()
    grid_points = np.stack((grid_points_x, grid_points_y), axis=1)
    for time_point in xrange(delays.size):
        print 'resample', time_point
        grid_values = data[time_point].ravel()
        interpolator = struc_3d.get_interpolator(grid_points, grid_values)
        values = interpolator(points[:, 0], points[:, 1])
        values_by_time_point.append(values)
    return points, values_by_time_point


def recombine_planes(planes, points_by_plane, values_by_plane, delays):
    r"""
    Recombine planes so that ug.save_to_unstructured_grid
    function can be used.
    """
    measurements = {'delays': delays,
                    'x_out': points_by_plane[0][:, 0],
                    'y_out': points_by_plane[0][:, 1],
                    'z_out': np.ones(points_by_plane[0][:, 0].size)*planes[0],
                    'a_out': values_by_plane[0]}
    for i, plane in enumerate(planes[1:]):
        measurements['x_out'].append(points_by_plane[i][:, 0])
        measurements['y_out'].append(points_by_plane[i][:, 1])
        measurements['z_out'].append(np.ones(points_by_plane[i].shape[0])*plane)
        measurements['a_out'].append(values_by_plane[i])
    return measurements


def parse_args():
    r"""
    """
    parser = argparse.ArgumentParser(description='Create unstructured VTK from measurements')
    parser.add_argument('--bx_extent',
                        help='spatial extent of Bx measurements',
                        nargs=6, type=float,
                        default=[-0.032, 0.028, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--bx_filter_extent',
                        help="spatial extent of interpolated grid"
                             "on which to filter Bx measurements",
                        nargs=4, type=float,
                        default=[-0.026, 0.025, -0.019, 0.029])
    parser.add_argument('--by_extent',
                        help='spatial extent of By measurements',
                        nargs=6, type=float,
                        default=[-0.032, 0.028, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--bxby_only',
                        help='flag to filter Bx and By only.',
                        default=False,
                        action='store_true')
    parser.add_argument('--bxby_extent',
                        help='spatial extent of Bx and By measurements',
                        nargs=6, type=float,
                        default=[-0.032, 0.026, -0.06, 0.043, 0.249, 0.416])
    parser.add_argument('--bxby_filter_extent',
                        help="spatial extent of interpolated grid"
                        "on which to filter Bx and By",
                        nargs=4, type=float,
                        default=[-0.032, 0.026, -0.06, 0.043, 0.249, 0.416])
    parser.add_argument('--bz_extent',
                        help='spatial extent of Bz measurements',
                        nargs=6, type=float,
                        default=[-0.032, 0.028, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--te_extent',
                        help='spatial extent of temperature measurements',
                        nargs=6, type=float,
                        default=[-0.026, 0.028, -0.03, 0.028, 0.249, 0.416])
    parser.add_argument('--te_bounds',
                        help='sensible bounds for temperature measurements',
                        nargs=2, type=float,
                        default=[1e-3, 1e3])
    parser.add_argument('--n_extent',
                        help='spatial extent of density measurements',
                        nargs=6, type=float,
                        default=[-0.026, 0.028, -0.03, 0.028, 0.249, 0.416])
    parser.add_argument('--n_bounds',
                        help='sensible bounds for density measurements',
                        nargs=2, type=float,
                        default=[1e3, 1e22])
    parser.add_argument('--mach_time_steps',
                        help='# of time steps to extract from one gyration', type=int,
                        default=250)
    parser.add_argument('--shot_database', help='path to shot database',
                        default='/home/jensv/rsx/jens_analysis/helicity_tools/shots_database/shots.db')
    parser.add_argument('--table_name', help='name of sql table',
                        default='Shots')
    parser.add_argument('--min_spectral',
                        help=("minimum spectral energy around gyration"
                              "frequency to include shot"),
                        type=float,
                        default=1.6e-8)
    parser.add_argument('--mach_y_extent',
                        help='spatial extent of mach measurements to include',
                        nargs=6, type=float,
                        default=[-0.052, 0.052, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--mach_z_extent',
                        help='spatial extent of mach measurements to include',
                        nargs=6, type=float,
                        default=[-0.032, 0.032, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--mach_bounds',
                        help='bounds on mach measurements', nargs=2, type=float,
                        default=[-10, 10])
    parser.add_argument('--output_prefix',
                        help='prefix of output files',
                        default='_filtered_unstructured_')
    parser.add_argument('--filter_spatial_increment',
                        help='spatial increment of interpolated grid for filtering',
                        default=0.0005, type=float)
    parser.add_argument('--filter_sigma',
                        help='standard deviation of gaussian filter',
                        type=float,
                        default=3)
    parser.add_argument('--filter_truncate',
                        help='truncate Gaussian filter at this multiple of sigma',
                        type=float,
                        default=3)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
