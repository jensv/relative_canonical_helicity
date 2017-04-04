#! /Users/vonderlinden2/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 18:07:14 2016

@author: Jens von der Linden
"""
import argparse
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from datetime import date
from datetime import datetime
import os

from write_to_vtk.read_unstructured_vtk import read_unstructured_vtk
from write_to_vtk import structured_3d_vtk as struc_3d

def main(args):
    r"""
    """
    just_magnetic = args.just_magnetic
    now = datetime.now().strftime("%Y-%m-%d-%H-%M")
    out_dir = '../output/' + args.output_prefix + '/' + now + '/'
    try:
        os.makedirs(out_dir)
    except:
        pass

    in_dir = args.input_path + args.input_date + '/'
    in_file = args.input_file_text

    for time_point in xrange(args.time_steps):
        print time_point
        time_str = str(time_point).zfill(4)
        bx_points, bx_values = read_unstructured_vtk(in_dir + 'bx' +
                                                     in_file + time_str + '.vtk')
        by_points, by_values = read_unstructured_vtk(in_dir + 'by' +
                                                     in_file + time_str + '.vtk')
        bz_points, bz_values = read_unstructured_vtk(in_dir + 'bz' +
                                                     in_file + time_str + '.vtk')
        if not just_magnetic:
            n_points, n_values = read_unstructured_vtk(in_dir + 'n' +
                                                       in_file + time_str + '.vtk')
            te_points, te_values = read_unstructured_vtk(in_dir + 'te' +
                                                         in_file + time_str + '.vtk')
            mach_y_points, mach_y_values = read_unstructured_vtk(in_dir + 'mach_y' +
                                                                 in_file + time_str + '.vtk')
            mach_z_points, mach_z_values = read_unstructured_vtk(in_dir + 'mach_z' +
                                                                 in_file + time_str + '.vtk')

        bx_interpolator = struc_3d.get_interpolator(bx_points, bx_values)
        by_interpolator = struc_3d.get_interpolator(by_points, by_values)
        bz_interpolator = struc_3d.get_interpolator(bz_points, bz_values)
        if not just_magnetic:
            te_interpolator = struc_3d.get_interpolator(te_points, te_values)
            n_interpolator = struc_3d.get_interpolator(n_points, n_values)
            mach_y_interpolator = struc_3d.get_interpolator(mach_y_points[:, :2], mach_y_values)
            mach_z_interpolator = struc_3d.get_interpolator(mach_z_points[:, :2], mach_z_values)

        (x_min, x_max,
         y_min, y_max,
         z_min, z_max) = args.joint_extent

        mesh = np.meshgrid(np.linspace(x_min, x_max,
                                       np.ceil((x_max-x_min)/
                                               args.spatial_increment)),
                           np.linspace(y_min, y_max,
                                       np.ceil((y_max-y_min)/
                                               args.spatial_increment)),
                           np.linspace(z_min, z_max,
                                       np.ceil((z_max-z_min)/
                                               args.spatial_increment)))


        bx_grad = struc_3d.triangulate_grad(mesh, bx_interpolator,
                                            increment=args.derivative_increment)
        by_grad = struc_3d.triangulate_grad(mesh, by_interpolator,
                                            increment=args.derivative_increment)
        bz_grad = struc_3d.triangulate_grad(mesh, bz_interpolator,
                                            increment=args.derivative_increment)
        if not just_magnetic:
            te_grad = struc_3d.triangulate_grad(mesh, te_interpolator,
                                                increment=args.derivative_increment)
            n_grad = struc_3d.triangulate_grad(mesh, n_interpolator,
                                               increment=args.derivative_increment)
            plane_mesh = [mesh[0][:, :, 0], mesh[1][:, :, 0]]
            mach_y_grad_plane = struc_3d.triangulate_grad(plane_mesh, mach_y_interpolator,
                                                          increment=args.derivative_increment)
            mach_z_grad_plane = struc_3d.triangulate_grad(plane_mesh, mach_z_interpolator,
                                                          increment=args.derivative_increment)


        bx, by, bz = struc_3d.vector_on_mesh((bx_interpolator,
                                              by_interpolator,
                                              bz_interpolator), mesh)
        bx, by, bz = struc_3d.add_vacuum_field([bx, by, bz],
                                               vacuum_field=args.bias_field_magnitude)
        if not just_magnetic:
            te = struc_3d.scalar_on_mesh(te_interpolator, mesh)
            n = struc_3d.scalar_on_mesh(n_interpolator, mesh)
            mach_y_plane = struc_3d.scalar_on_mesh(mach_y_interpolator,
                                               plane_mesh)
            mach_z_plane= struc_3d.scalar_on_mesh(mach_z_interpolator,
                                              plane_mesh)

            mach_y = np.repeat(mach_y_plane[:, :, np.newaxis],
                           mesh[0].shape[2], axis=2)
            mach_z = np.repeat(mach_z_plane[:, :, np.newaxis],
                           mesh[0].shape[2], axis=2)

            mach_y_dx = np.repeat(mach_y_grad_plane[0][:, :, np.newaxis],
                              mesh[0].shape[2], axis=2)
            mach_y_dy = np.repeat(mach_y_grad_plane[1][:, :, np.newaxis],
                              mesh[0].shape[2], axis=2)
            mach_y_dz = np.zeros(mesh[0].shape)

            mach_z_dx = np.repeat(mach_z_grad_plane[0][:, :, np.newaxis],
                              mesh[0].shape[2], axis=2)
            mach_z_dy = np.repeat(mach_z_grad_plane[0][:, :, np.newaxis],
                              mesh[0].shape[2], axis=2)
            mach_z_dz = np.zeros(mesh[0].shape)

            fields = ([bx] + [by] + [bz] + [n] + [te] +
                      [mach_y] + [mach_z] +
                      list(bx_grad) +
                      list(by_grad) +
                      list(bz_grad) +
                      list(n_grad) +
                      list(te_grad) +
                      [mach_y_dx] + [mach_y_dy] + [mach_y_dz] +
                      [mach_z_dx] + [mach_z_dy] + [mach_z_dz])

            quantity_names = ['B_x', 'B_y', 'B_z',
                              'n', 'Te',
                              'mach_y', 'mach_z',
                              'B_x_dx', 'B_x_dy', 'B_x_dz',
                              'B_y_dx', 'B_y_dy', 'B_y_dz',
                              'B_z_dx', 'B_z_dy', 'B_z_dz',
                              'n_dx', 'n_dy', 'n_dz',
                              'Te_dx', 'Te_dy', 'Te_dz',
                              'mach_y_dx', 'mach_y_dy', 'mach_y_dz',
                              'mach_z_dx', 'mach_z_dy', 'mach_z_dz']
        else:
            fields = ([bx] + [by] + [bz] +
                      list(bx_grad) +
                      list(by_grad) +
                      list(bz_grad))

            quantity_names = ['B_x', 'B_y', 'B_z',
                              'B_x_dx', 'B_x_dy', 'B_x_dz',
                              'B_y_dx', 'B_y_dy', 'B_y_dz',
                              'B_z_dx', 'B_z_dy', 'B_z_dz']


        x, y, z, variables = struc_3d.prepare_for_rectilinear_grid(mesh, fields,
                                                                   quantity_names)

        vtk_file_path = out_dir + args.output_prefix
        struc_3d.write_fields_to_rectilinear_grid(vtk_file_path,
                                                  x, y, z, variables,
                                                  time_point)

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
                        default='2017-04-04-13-44')
    parser.add_argument('--input_file_text',
                        help='input file name',
                        default='_boxed_unstructured_')
    parser.add_argument('--spatial_increment',
                        help='Spatial increment of output file grids',
                        type=float, default=0.001)
    parser.add_argument('--derivative_increment',
                        help=("spatial increment used to determine"
                              "tetrahedron derivative of Delaunay"),
                        type=float, default=0.0001)
    parser.add_argument('--joint_extent',
                        help='overlapping spatial extent of all parameters',
                        nargs=6, type=float,
                        default=[-0.022, 0.024, -0.02, 0.018, 0.249, 0.416])
    parser.add_argument('--output_prefix',
                        help='prefix of output files',
                        default='data_interp_to_rect_grid')
    parser.add_argument('--bias_field_magnitude',
                        help='magnitude of axial bias magnetic field',
                        type=float,
                        default=0.02)
    parser.add_argument('--time_steps',
                        help='number of time steps', type=int,
                        default=250)
    parser.add_argument('--just_magnetic',
                        help='only interpolate bdot measurements',
                        action='store_true', default=False)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
