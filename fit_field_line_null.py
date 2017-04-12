r"""
Fits field line null.

Created March 28 2017 by Jens von der Linden.
"""

import argparse
from scipy.interpolate import LinearNDInterpolator
import numpy as np
from scipy.optimize import leastsq
from scipy import odr
from scipy import ndimage
from scipy.integrate import odeint, dblquad
from datetime import datetime
import os
from scipy.interpolate import LinearNDInterpolator

from write_to_vtk.read_unstructured_vtk import read_unstructured_vtk
from write_to_vtk import structured_3d_vtk as struc_3d


def main(args):
    r"""
    """
    now = datetime.now().strftime("%Y-%m-%d-%H-%M")
    out_dir = '../output/' + args.output_prefix + '/' + now + '/'
    try:
        os.makedirs(out_dir)
    except:
        pass

    in_dir = args.input_path + args.input_date + '/'
    in_file = args.input_file_text

    centroids = []


    bxby_extents = {0: args.bxby_extent_0,
                    1: args.bxby_extent_1,
                    2: args.bxby_extent_2,
                    3: args.bxby_extent_3}
    bz_extents = {0: args.bz_extent_0,
                  1: args.bz_extent_1,
                  2: args.bz_extent_2,
                  3: args.bz_extent_3}

    bxby_extent = bxby_extents[args.plane_number]
    bz_extent = bz_extents[args.plane_number]

    for time_point in xrange(args.time_steps):
        print time_point
        time_str = str(time_point).zfill(4)
        bx_points, bx_values = read_unstructured_vtk(in_dir + 'bx' +
                                                     in_file + time_str + '.vtk')
        by_points, by_values = read_unstructured_vtk(in_dir + 'by' +
                                                     in_file + time_str + '.vtk')
        bz_points, bz_values = read_unstructured_vtk(in_dir + 'bz' +
                                                     in_file + time_str + '.vtk')
        z_value = np.unique(bx_points[:, 2])[args.plane_number]

        z_index = np.where(bx_points[:, 2] == z_value)[0]
        bx_points = bx_points[z_index, :-1]
        bx_values = bx_values[z_index]
        z_index = np.where(by_points[:, 2] == z_value)[0]
        by_points = by_points[z_index, :-1]
        by_values = by_values[z_index]
        z_index = np.where(bz_points[:, 2] == z_value)[0]
        bz_points = bz_points[z_index, :-1]
        bz_values = bz_values[z_index]

        bx_interpolator = struc_3d.get_interpolator(bx_points, bx_values)
        by_interpolator = struc_3d.get_interpolator(by_points, by_values)
        bz_interpolator = struc_3d.get_interpolator(bz_points, bz_values)
        grid_extent = [bxby_extent[0], bxby_extent[1],
                       -0.02, bxby_extent[3]]
        grid = np.meshgrid(np.linspace(grid_extent[0], grid_extent[1],
                                       (grid_extent[1] - grid_extent[0])/
                                       args.spatial_increment),
                           np.linspace(grid_extent[2], grid_extent[3],
                                       (grid_extent[3] - grid_extent[2])/
                                       args.spatial_increment))

        (centroid, center_points,
         radii, streamlines,
         max_index) = find_field_null(grid,
                                      bx_interpolator,
                                      by_interpolator,
                                      launch_point_step_factor=0.05,
                                      integration_length=20)
        centroids.append(centroid)

    centroids = np.asarray(centroids)
    np.savetxt(out_dir + '/field_nulls.txt', centroids,
               header=("magnetic field null positions in z plane # %d plane,"
                       "determined by"
                       "fitting circles to integrated field lines starting at max"
                       "magnitude and moving successive towards the center of circles."
                       % args.plane_number))


def d_l(l, t, interpolator_x, interpolator_y):
    r"""
    Returns d_l for the field line integrator.
    """
    return np.asarray([interpolator_x([l[0], l[1]])[0],
                       interpolator_y([l[0], l[1]])[0]])


def to_min(params, points):
    r"""
    Returns circle expression to minimize with least squares.
    """
    a = 2.*params[0]
    b = 2.*params[1]
    c = params[2]**2 - params[1]**2 - params[0]**2
    return a*points[0] + b*points[1] + c - points[0]**2 - points[1]**2


def find_field_null(grid, bx_interpolator, by_interpolator,
                    distance_thres=0.001, filter_size=2,
                    integration_length=10, integration_steps=100,
                    launch_point_step_factor=0.1, max_count=50,
                    params_guess=[0, 0, 0.01]):
    r"""
    Find Bx-By field null in a x-y plane
    by integrating field lines and fitting a circle to them.
    Move towards the center and iterate process.
    If leaving the measurement plane
    extrapolate from last fit circle.
    Start close to the Bx-By field max.
    """
    b_fields_x = bx_interpolator(grid[0][:, :], grid[1][:, :])
    b_fields_y = by_interpolator(grid[0][:, :], grid[1][:, :])
    b_fields = [b_fields_x, b_fields_y]
    x_min, x_max = grid[0].min(), grid[0].max()
    y_min, y_max = grid[1].min(), grid[1].max()
    magnitude = np.sqrt(b_fields[0][:, :]**2 + b_fields[1][:, :]**2)
    filtered_magnitude = ndimage.gaussian_filter(magnitude, filter_size)
    max_index = np.unravel_index( np.nanargmax(filtered_magnitude),
                                 filtered_magnitude.shape)
    center_points = []
    radii = []
    center_points = []
    streamlines = []
    direction = [0, 0]
    distance = 100
    launch_point = (grid[0][:][max_index], grid[1][:][max_index])
    count = 0
    while distance >= distance_thres:
        #print 'launch', launch_point
        #print distance
        t2 = np.linspace(0, integration_length, integration_steps)
        t1 = np.linspace(0, -integration_length, integration_steps)
        stream2 = odeint(d_l, launch_point, t2, args=(bx_interpolator, by_interpolator))
        stream1 = odeint(d_l, launch_point, t1, args=(bx_interpolator, by_interpolator))
        streamline = np.concatenate((stream1, stream2))
        size0 = np.sum(np.invert(np.isnan(streamline[:, 0])))
        size1 = np.sum(np.invert(np.isnan(streamline[:, 1])))
        min_index = np.argmin([size0, size1])
        min_size = [size0, size1][min_index]
        streamline = streamline[np.invert(np.isnan(streamline[:, min_index]))].reshape(min_size, 2)
        circle_params, success = leastsq(to_min, params_guess,
                                         args=np.asarray([streamline[:, 0],
                                                          streamline[:, 1]]))
        direction = [circle_params[0] - launch_point[0], circle_params[1] - launch_point[1]]
        distance = np.sqrt(direction[0]**2. + direction[1]**2.)
        center_point = (circle_params[0], circle_params[1])
        launch_point = [launch_point[0] + direction[0] * launch_point_step_factor,
                        launch_point[1] + direction[1] * launch_point_step_factor]
        center_points.append(center_point)
        #print 'center', center_point
        radii.append(circle_params[0])
        streamlines.append(streamline)
        if (launch_point[0] <= x_min or
            launch_point[0] >= x_max or
            launch_point[1] <= y_min or
            launch_point[1] >= y_max or
            count > max_count):
            break
        count += 1
    field_null = center_point
    return field_null, center_points, radii, streamlines, max_index


def integrate_flux(centroid, radius, bz_interpolator, limits, bias_field=0.02):
    r"""
    Return axial magnetic flux and error estimate integrated
    in a circle of given radius around a given centroid.
    """
    if (centroid[0] - radius < limits[0] or centroid[0] + radius > limits[1] or
        centroid[1] - radius < limits[2] or centroid[1] + radius > limits[3]):
        return -1
    gfun = lambda x: -np.sqrt(radius**2 - (x-centroid[0])**2)
    hfun = lambda x: np.sqrt(radius**2 - (x-centroid[0])**2)
    bz_interpolator_bias  = lambda x, y: bz_interpolator(x, y) + bias_field
    return dblquad(bz_interpolator_bias, centroid[0] - radius,
                   centroid[0] + radius, gfun, hfun)


def parse_args():
    r"""
    """
    parser = argparse.ArgumentParser(description=("Create VTK files of"
                                                  "interpolated measurements"))
    parser.add_argument('--input_path',
                        help='path to input files',
                        default='../output/boxed_unstructured_measurements/')
    parser.add_argument('--input_date',
                        help='time stamp of input files',
                        default='2017-04-11-21-07')
    parser.add_argument('--input_file_text',
                        help='input file name',
                        default='_boxed_unstructured_')
    parser.add_argument('--spatial_increment',
                        help='Spatial increment of grids',
                        type=float, default=0.001)
    parser.add_argument('--derivative_increment',
                        help=("spatial increment used to determine"
                              "tetrahedron derivative of Delaunay"),
                        type=float, default=0.0001)
    parser.add_argument('--plane_number',
                        help="z-Plane number on which to find field nulls.",
                        type=int, default=0)
    parser.add_argument('--bxby_extent_0',
                        help='overlapping spatial extent of bx by',
                        nargs=6, type=float,
                        default=[-0.027, 0.025, -0.057, 0.041, 0.249, 0.416])
    parser.add_argument('--bxby_extent_1',
                        help='overlapping spatial extent of bx by',
                        nargs=6, type=float,
                        default=[-0.027, 0.027, -0.073, 0.041, 0.249, 0.416])
    parser.add_argument('--bxby_extent_2',
                        help='overlapping spatial extent of bx by',
                        nargs=6, type=float,
                        default=[-0.047, 0.031, -0.021, 0.028, 0.249, 0.416])
    parser.add_argument('--bxby_extent_3',
                        help='overlapping spatial extent of bx by',
                        nargs=6, type=float,
                        default=[-0.061, 0.031, -0.026, 0.03, 0.249, 0.416])
    parser.add_argument('--bz_extent_0',
                        help='spatial extent of bz',
                        nargs=6, type=float,
                        default=[-0.027, 0.025, -0.06, 0.041, 0.249, 0.416])
    parser.add_argument('--bz_extent_1',
                        help='spatial extent of bz',
                        nargs=6, type=float,
                        default=[-0.27, 0.027, -0.076, 0.041, 0.249, 0.416])
    parser.add_argument('--bz_extent_2',
                        help='spatial extent of bz',
                        nargs=6, type=float,
                        default=[-0.044, 0.031, -0.021, 0.03, 0.249, 0.416])
    parser.add_argument('--bz_extent_3',
                        help='spatial extent of bz',
                        nargs=6, type=float,
                        default=[-0.072, 0.031, -0.026, 0.03, 0.249, 0.416])
    parser.add_argument('--output_prefix',
                        help='prefix of output files',
                        default='field_nulls')
    parser.add_argument('--bias_field_magnitude',
                        help='magnitude of axial bias magnetic field',
                        type=float,
                        default=0.02)
    parser.add_argument('--time_steps',
                        help='number of time steps', type=int,
                        default=250)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    main(args)
