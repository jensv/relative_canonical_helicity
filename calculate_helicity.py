#! /Users/vonderlinden2/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 18:07:14 2016

@author: Jens von der Linden
"""
import argparse
import numpy as np
from scipy.constants import mu_0
from scipy.constants import elementary_charge as q_e
from scipy.constants import proton_mass as m_i
from scipy.integrate import trapz
from datetime import date
from datetime import datetime
import visit_writer
import os
from vector_comparison import read_rectilinear_vtk as rrv
from write_to_vtk import structured_3d_vtk as struc_3d
import write_to_vtk.helicities as hel


def main(args):
    r"""
    """
    just_magnetic = args.just_magnetic
    now = datetime.now().strftime("%Y-%m-%d-%H-%M")
    out_dir = '../output/' + args.output_prefix + '/' + now
    try:
        os.makedirs(out_dir)
    except:
        pass
    in_dir = args.input_path + args.input_date + '/'
    in_file = args.input_file_text
    input_path = in_dir + in_file
    xlow, xhigh = args.to_cut[:2]
    ylow, yhigh = args.to_cut[2:4]
    zlow, zhigh = args.to_cut[4:6]

    if not just_magnetic:
        fields_to_read = ('B_x', 'B_y', 'B_z',
                          'B_ref_x', 'B_ref_y', 'B_ref_z',
                          'A_x', 'A_y', 'A_z',
                          'A_ref_x', 'A_ref_y', 'A_ref_z',
                          'u_i_x_plus', 'u_i_y', 'u_i_z',
                          'u_i_ref_x', 'u_i_ref_y', 'u_i_ref_z',
                          'w_i_raw_x_plus', 'w_i_raw_y_plus', 'w_i_raw_z_plus',
                          'w_i_raw_ref_x', 'w_i_raw_ref_y', 'w_i_raw_ref_z',
                          'w_i_x_plus', 'w_i_y_plus', 'w_i_z_plus',
                          'w_i_ref_x', 'w_i_ref_y', 'w_i_ref_z', 'n')
    else:
        fields_to_read = ('B_x', 'B_y', 'B_z',
                          'B_ref_x', 'B_ref_y', 'B_ref_z',
                          'A_x', 'A_y', 'A_z',
                          'A_ref_x', 'A_ref_y', 'A_ref_z')


    for time_point in xrange(args.time_steps):
        print time_point
        input_file = (input_path +
                      str(time_point).zfill(4) +
                      '.vtk')
        fields = {}
        for field in fields_to_read:
            read = rrv.read_scalar(input_file, field)
            fields[field] = read[1][ylow:yhigh, xlow:xhigh, zlow:zhigh]
        mesh = read[0]
        mesh[0] = mesh[0][ylow:yhigh, xlow:xhigh, zlow:zhigh]
        mesh[1] = mesh[1][ylow:yhigh, xlow:xhigh, zlow:zhigh]
        mesh[2] = mesh[2][ylow:yhigh, xlow:xhigh, zlow:zhigh]
 
        dx = mesh[0][0, 1, 0] - mesh[0][0, 0, 0]
        dy = mesh[1][1, 0, 0] - mesh[1][0, 0, 0]
        dz = mesh[2][0, 0, 1] - mesh[2][0, 0, 0]

        b_field = (fields['B_x'],
                   fields['B_y'],
                   fields['B_z'])
        b_field_ref = (fields['B_ref_x'],
                       fields['B_ref_y'],
                       fields['B_ref_z'])
        vector_potential = (fields['A_x'],
                            fields['A_y'],
                            fields['A_z'])
        vector_potential_ref = (fields['A_ref_x'],
                                fields['A_ref_y'],
                                fields['A_ref_z'])

        if not just_magnetic:
            ion_velocity = (fields['u_i_x_plus'],
                            fields['u_i_y'],
                            fields['u_i_z'])
            ion_velocity_ref = (fields['u_i_ref_x'],
                                fields['u_i_ref_y'],
                                fields['u_i_ref_z'])
            ion_vorticity = (fields['w_i_x_plus'],
                             fields['w_i_y_plus'],
                             fields['w_i_z_plus'])
            ion_vorticity_raw = (fields['w_i_raw_x_plus'],
                                 fields['w_i_raw_y_plus'],
                                 fields['w_i_raw_z_plus'])
            ion_vorticity_ref = (fields['w_i_ref_x'],
                                 fields['w_i_ref_y'],
                                 fields['w_i_ref_z'])
            ion_vorticity_raw_ref =(fields['w_i_raw_ref_x'],
                                    fields['w_i_raw_ref_y'],
                                    fields['w_i_raw_ref_z'])
            n = fields['n']

            kin_helicity = hel.kinetic_helicity(ion_velocity,
                                                ion_vorticity,
                                                dx, dy, dz)
            kin_helicity_raw_vorticity = hel.kinetic_helicity(ion_velocity,
                                                              ion_vorticity_raw,
                                                              dx, dy, dz)
            cross_helicity = hel.cross_helicity(ion_velocity,
                                                b_field,
                                                dx, dy, dz)
            rel_kin_helicity = hel.rel_kin_helicity(ion_velocity,
                                                    ion_velocity_ref,
                                                    ion_vorticity,
                                                    ion_vorticity_ref,
                                                    dx, dy, dz)
            rel_cross_helicity = hel.rel_cross_helicity(ion_velocity,
                                                        ion_velocity_ref,
                                                        b_field,
                                                        b_field_ref,
                                                        dx, dy, dz)
            rel_kin_helicity_raw_vorticity = hel.rel_kin_helicity(ion_velocity,
                                                                  ion_velocity_ref,
                                                                  ion_vorticity_raw,
                                                                  ion_vorticity_raw_ref,
                                                                  dx, dy, dz)
            rel_cross_helicity_raw_vorticity = hel.rel_cross_helicity(ion_velocity,
                                                                      ion_velocity_ref,
                                                                      b_field,
                                                                      b_field_ref,
                                                                      dx, dy, dz)

            ## helicities with density dependence
            ##
            kin_helicity_n = hel.kinetic_helicity(ion_velocity,
                                                  ion_vorticity,
                                                  dx, dy, dz, density=n)
            kin_helicity_raw_vorticity_n = hel.kinetic_helicity(ion_velocity,
                                                                ion_vorticity_raw,
                                                                dx, dy, dz,
                                                                density=n)
            cross_helicity_n = hel.cross_helicity(ion_velocity,
                                                  b_field,
                                                  dx, dy, dz,
                                                  density=n)
            rel_kin_helicity_n = hel.rel_kin_helicity(ion_velocity,
                                                      ion_velocity_ref,
                                                      ion_vorticity,
                                                      ion_vorticity_ref,
                                                      dx, dy, dz, density=n)
            rel_cross_helicity_n = hel.rel_cross_helicity(ion_velocity,
                                                          ion_velocity_ref,
                                                          b_field,
                                                          b_field_ref,
                                                          dx, dy, dz, density=n)
            rel_kin_helicity_raw_vorticity_n = hel.rel_kin_helicity(ion_velocity,
                                                                    ion_velocity_ref,
                                                                    ion_vorticity_raw,
                                                                    ion_vorticity_raw_ref,
                                                                    dx, dy, dz, density=n)
            rel_cross_helicity_raw_vorticity_n = hel.rel_cross_helicity(ion_velocity,
                                                                        ion_velocity_ref,
                                                                        b_field,
                                                                        b_field_ref,
                                                                        dx, dy, dz,
                                                                        density=n)
            mag_helicity_n = hel.magnetic_helicity(vector_potential,
                                                   b_field,
                                                   dx, dy, dz,
                                                   density=n)

            rel_mag_helicity_n = hel.rel_mag_helicity(vector_potential,
                                                      vector_potential_ref,
                                                      b_field,
                                                      b_field_ref,
                                                      dx, dy, dz,
                                                      density=n)



        mag_helicity = hel.magnetic_helicity(vector_potential,
                                             b_field,
                                             dx, dy, dz)

        rel_mag_helicity = hel.rel_mag_helicity(vector_potential,
                                                vector_potential_ref,
                                                b_field,
                                                b_field_ref,
                                                dx, dy, dz)
        if not just_magnetic:
            with open(out_dir +
                      "/kinetic_helicity.txt", "a") as myfile:
                myfile.write(str(kin_helicity) + '\n')
            with open(out_dir +
                      "/kinetic_helicity_raw_vorticity.txt",
                      "a") as myfile:
                myfile.write(str(kin_helicity_raw_vorticity) + '\n')
            with open(out_dir +
                      "/cross_helicity.txt", "a") as myfile:
                myfile.write(str(cross_helicity) + '\n')
            with open(out_dir +
                      "/relative_kinetic_helicity.txt", "a") as myfile:
                myfile.write(str(rel_kin_helicity) + '\n')
            with open(out_dir +
                      "/relative_cross_helicity.txt", "a") as myfile:
                myfile.write(str(rel_cross_helicity) + '\n')
            with open(out_dir +
                      "/relative_kinetic_helicity_raw_vorticity.txt",
                      "a") as myfile:
                myfile.write(str(rel_kin_helicity_raw_vorticity) + '\n')
            with open(out_dir +
                      "/relative_cross_helicity_raw_vorticity.txt",
                      "a") as myfile:
                myfile.write(str(rel_cross_helicity_raw_vorticity) + '\n')

           ## with density dependence
           ##
            with open(out_dir +
                      "/kinetic_helicity_n_dependence.txt", "a") as myfile:
                myfile.write(str(kin_helicity_n) + '\n')
            with open(out_dir +
                      "/kinetic_helicity_raw_vorticity_n_dependence.txt",
                      "a") as myfile:
                myfile.write(str(kin_helicity_raw_vorticity_n) + '\n')
            with open(out_dir +
                      "/cross_helicity_n_dependence.txt", "a") as myfile:
                myfile.write(str(cross_helicity_n) + '\n')
            with open(out_dir +
                      "/relative_kinetic_helicity_n_dependence.txt", "a") as myfile:
                myfile.write(str(rel_kin_helicity_n) + '\n')
            with open(out_dir +
                      "/relative_cross_helicity_n_dependence.txt", "a") as myfile:
                myfile.write(str(rel_cross_helicity_n) + '\n')
            with open(out_dir +
                      "/relative_kinetic_helicity_raw_vorticity_n_dependence.txt",
                      "a") as myfile:
                myfile.write(str(rel_kin_helicity_raw_vorticity_n) + '\n')
            with open(out_dir +
                      "/relative_cross_helicity_raw_vorticity_n_dependence.txt",
                      "a") as myfile:
                myfile.write(str(rel_cross_helicity_raw_vorticity_n) + '\n')
            with open(out_dir +
                      "/magnetic_helicity_n_dependence.txt", "a") as myfile:
                myfile.write(str(mag_helicity_n) + '\n')
            with open(out_dir +
                      "/relative_magnetic_helicity_n_dependence.txt", "a") as myfile:
                myfile.write(str(rel_mag_helicity_n) + '\n')


        with open(out_dir +
                  "/magnetic_helicity.txt", "a") as myfile:
            myfile.write(str(mag_helicity) + '\n')
        with open(out_dir +
                  "/relative_magnetic_helicity.txt", "a") as myfile:
            myfile.write(str(rel_mag_helicity) + '\n')

def parse_args():
    r"""
    """
    parser = argparse.ArgumentParser(description='Calculate helicity')
    parser.add_argument('--input_path',
                        help='path to input files',
                        default='../output/canonical_quantities/')
    parser.add_argument('--input_date',
                        help='time stamp of input files',
                        default='2017-04-04-16-33')
    parser.add_argument('--input_file_text',
                        help='input file name',
                        default='canonical_quantities')
    parser.add_argument('--output_prefix',
                        help='prefix of output files',
                        default='helicity')
    parser.add_argument('--time_steps',
                        help='number of time steps',
                        type=int,
                        default=250)
    parser.add_argument('--to_cut',
                        help=("number of indicies to cut on"
                              "each side to remove NaNs e.g."
                              "xlow xhigh ylow yhigh zlow zhigh"),
                        nargs=6,
                        default=['0', 'None', '0', 'None', '0', 'None'])
    parser.add_argument('--just_magnetic',
                        help='only calculate magnetic quantities',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    to_cut_list = []
    for to_cut in args.to_cut:
        if to_cut == 'None':
            to_cut_list.append(None)
        else:
            to_cut_list.append(int(to_cut))
    args.to_cut = to_cut_list
    return args


if __name__ == '__main__':
    args = parse_args()
    main(args)
