#! /Users/vonderlinden2/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on Feb 14 2017

@author: Jens von der Linden

Read processed measurements
and write to unstructured vtk files files.
The processed measurements are magnetic field, density, temperature
and Mach number.
Mach number is calculated directly from the raw MDSplus Mach probe
voltages, all other data is read from IDL output files from the RSX data anlysis
scripts.
Only measurement points slightly larger than the
intended interpolation spaces are kept and Delaunay triangulated
to define the unstructured grid.
"""
import argparse
import numpy as np
from datetime import date
from datetime import datetime
import os

from mach_probe_analysis import ion_current_to_mach_number as ic_to_mach
from read_from_sql import read_from_sql
from write_to_vtk import prepare_measurements as pm
from write_to_vtk import unstructured_grid as ug


def main(args):
    r"""
    Read processed measurements and write to unstructured vtk files.
    """
    now = datetime.now().strftime("%Y-%m-%d-%H-%M")
    out_dir = '../output/boxed_unstructured_measurements/' + now
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

    ug.save_to_unstructured_grid(bx_all_planes, 'bx', out_dir,
                                 prefix=args.output_prefix)
    ug.save_to_unstructured_grid(by_all_planes, 'by', out_dir,
                                 prefix=args.output_prefix)
    ug.save_to_unstructured_grid(bz_all_planes, 'bz', out_dir,
                                 prefix=args.output_prefix)
    ug.save_to_unstructured_grid(te_three_planes, 'te', out_dir,
                                 prefix=args.output_prefix)
    ug.save_to_unstructured_grid(n_three_planes, 'n', out_dir,
                                 prefix=args.output_prefix)
    ug.save_to_unstructured_grid(mach_y_plane, 'mach_y', out_dir,
                                 prefix=args.output_prefix)
    ug.save_to_unstructured_grid(mach_z_plane, 'mach_z', out_dir,
                                 prefix=args.output_prefix)

def parse_args():
    r"""
    Read Arguments.
    """
    parser = argparse.ArgumentParser(description='Create unstructured VTK from measurements')
    parser.add_argument('--bx_extent',
                        help='spatial extent of Bx measurements',
                        nargs=6, type=float,
                        default=[-0.032, 0.028, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--by_extent',
                        help='spatial extent of By measurements',
                        nargs=6, type=float,
                        default=[-0.032, 0.028, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--bz_extent',
                        help='spatial extent of Bz measurements',
                        nargs=6, type=float,
                        default=[-0.032, 0.028, -0.022, 0.032, 0.249, 0.416])
    parser.add_argument('--te_extent',
                        help='spatial extent of temperature measurements',
                        nargs=6, type=float,
                        default=[-0.026, 0.028, -0.03, 0.029, 0.249, 0.416])
    parser.add_argument('--te_bounds',
                        help='sensible bounds for temperature measurements',
                        nargs=2, type=float,
                        default=[1e-3, 1e3])
    parser.add_argument('--n_extent',
                        help='spatial extent of density measurements',
                        nargs=6, type=float,
                        default=[-0.026, 0.028, -0.03, 0.029, 0.249, 0.416])
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
                        default='_boxed_unstructured_')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    main(args)
