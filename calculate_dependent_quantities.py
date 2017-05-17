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
from astropy.convolution import convolve, convolve_fft
from scipy.ndimage import gaussian_filter
from scipy import ndimage
from scipy.signal import fftconvolve
import scipy
from datetime import date
from datetime import datetime
import visit_writer
import os

from vector_comparison import read_rectilinear_vtk as rrv
from write_to_vtk import structured_3d_vtk as struc_3d
from invert_curl.invert_curl import devore_invert_curl
from write_to_vtk.reference_field import determine_reference_fields as ref_fields
from vector_calculus import vector_calculus as vc

def main(args):
    r"""
    """
    just_magnetic = args.just_magnetic
    interpolate_nan = args.interpolate_nan
    now = datetime.now().strftime("%Y-%m-%d-%H-%M")
    out_dir = '../output/' + args.output_prefix + '/' + now + '/'
    try:
        os.makedirs(out_dir)
    except:
        pass
    in_dir = args.input_path + args.input_date + '/'
    in_file = args.input_file_text
    input_path = in_dir + in_file

    if not just_magnetic:
        fields_to_read = ('B_x', 'B_y', 'B_z',
                          'n', 'Te',
                          'mach_y', 'mach_z',
                          'B_x_dx', 'B_x_dy', 'B_x_dz',
                          'B_y_dx', 'B_y_dy', 'B_y_dz',
                          'B_z_dx', 'B_z_dy', 'B_z_dz',
                          'n_dx', 'n_dy', 'n_dz',
                          'Te_dx', 'Te_dy', 'Te_dz',
                          'mach_y_dx', 'mach_y_dy', 'mach_y_dz',
                          'mach_z_dx', 'mach_z_dy', 'mach_z_dz')
        mach_fields_to_read = ('mach_y', 'mach_y_dx', 'mach_y_dy', 'mach_y_dz')
    else:
        fields_to_read = ('B_x', 'B_y', 'B_z',
                          'B_x_dx', 'B_x_dy', 'B_x_dz',
                          'B_y_dx', 'B_y_dy', 'B_y_dz',
                          'B_z_dx', 'B_z_dy', 'B_z_dz')

    u_x_time_points_plus = np.roll(np.arange(250), int(np.round(250*0.25)))
    u_x_time_points_minus = np.roll(np.arange(250), -int(np.round(250*0.25)))
    xlow, xhigh, ylow, yhigh, zlow, zhigh = args.to_cut


    for time_point in xrange(args.time_steps):
        print time_point
        input_file = input_path + str(time_point).zfill(4) + '.vtk'
        fields = {}
        for field in fields_to_read:
            read = rrv.read_scalar(input_file, field)
            fields[field] = read[1][ylow:yhigh, xlow:xhigh, zlow:zhigh]
            if args.filter_gaussian:
                fields[field] = filter_data(fields[field],
                                            filter_sigma=args.filter_sigma,
                                            filter_truncate=args.filter_truncate)

        mesh = read[0]
        mesh[0] = mesh[0][ylow:yhigh, xlow:xhigh, zlow:zhigh]
        mesh[1] = mesh[1][ylow:yhigh, xlow:xhigh, zlow:zhigh]
        mesh[2] = mesh[2][ylow:yhigh, xlow:xhigh, zlow:zhigh]

        bx_grad = [fields['B_x_dx'], fields['B_x_dy'], fields['B_x_dz']]
        by_grad = [fields['B_y_dx'], fields['B_y_dy'], fields['B_y_dz']]
        bz_grad = [fields['B_z_dx'], fields['B_z_dy'], fields['B_z_dz']]

        ## current
        ##
        current = 1/mu_0*np.asarray(vc.curl(None,
                                            vector_grad=[bx_grad,
                                                         by_grad,
                                                         bz_grad]))
        current_smooth = []
        for direction in xrange(len(current)):
            current_dir = np.array(current[direction])
            for z_step in xrange(len(current[0, 0, :])):
                current_dir[:, :, z_step] = (gaussian_filter(current_dir[:, :, z_step],
                                                             args.filter_width,
                                                             mode='reflect'))
            current_smooth.append(current_dir)

        ## density and temperature
        ##
        if not just_magnetic:
            density = fields['n']
            grad_density = (fields['n_dx'], fields['n_dy'], fields['n_dz'])
            temperature = fields['Te']
            grad_temperature = (fields['Te_dx'], fields['Te_dy'], fields['Te_dz'])

            density_plane_normalized = normalize_scalar_by_plane(density)
            temperature_plane_normalized = normalize_scalar_by_plane(temperature)
            density_smooth = boxcar_filter(density,
                                           args.filter_width)
            temperature_smooth = boxcar_filter(temperature,
                                               args.filter_width)
            density_smooth_plane_normalized = normalize_scalar_by_plane(density_smooth)
            temperature_smooth_plane_normalized = normalize_scalar_by_plane(temperature_smooth)
            density_constant = args.density_constant_factor*np.ones(density.shape)
            temperature_smooth_norm = temperature_smooth/np.nanmax(temperature_smooth)
            temperature_norm = (temperature/
                                np.nanmax(temperature))

            ## velocity and vorticity
            ##
            mach_file = (input_path +
                         str(u_x_time_points_plus[time_point]).zfill(4) +
                         '.vtk')
            mach_x_fields = {}
            for field in mach_fields_to_read:
                read = rrv.read_scalar(mach_file, field)
                mach_x_fields[field] = read[1][ylow:yhigh, xlow:xhigh, zlow:zhigh]
            mach_x_p = mach_x_fields['mach_y']
            grad_mach_x_p = (mach_x_fields['mach_y_dx'],
                             mach_x_fields['mach_y_dy'],
                             mach_x_fields['mach_y_dz'])
            mach_file = (input_path +
                         str(u_x_time_points_minus[time_point]).zfill(4) +
                         '.vtk')
            mach_x_fields = {}
            for field in mach_fields_to_read:
                read = rrv.read_scalar(mach_file, field)
                mach_x_fields[field] = read[1][ylow:yhigh, xlow:xhigh, zlow:zhigh]
            mach_x_m = mach_x_fields['mach_y']
            grad_mach_x_m = (mach_x_fields['mach_y_dx'],
                             mach_x_fields['mach_y_dy'],
                             mach_x_fields['mach_y_dz'])
            mach_y = fields['mach_y']
            grad_mach_y = (fields['mach_y_dx'],
                           fields['mach_y_dy'],
                           fields['mach_y_dz'])
            mach_z = fields['mach_z']
            grad_mach_z = (fields['mach_z_dx'],
                           fields['mach_z_dy'],
                           fields['mach_z_dz'])
            te = temperature
            u_i_x_plus = np.sqrt(te*q_e/m_i)*mach_x_p
            u_i_x_minus = np.sqrt(te*q_e/m_i)*mach_x_m
            u_i_y = np.sqrt(te*q_e/m_i)*mach_y
            u_i_z = np.sqrt(te*q_e/m_i)*mach_z
            ion_velocity_p = [u_i_x_plus, u_i_y, u_i_z]
            ion_velocity_m = [u_i_x_minus, u_i_y, u_i_z]
            grad_ion_velocity_x_p = grad_ion_velocity(mach_x_p, grad_mach_x_p,
                                                      temperature, grad_temperature)
            grad_ion_velocity_x_m = grad_ion_velocity(mach_x_m, grad_mach_x_m,
                                                      temperature, grad_temperature)
            grad_ion_velocity_y = grad_ion_velocity(mach_y, grad_mach_y,
                                                    temperature, grad_temperature)
            grad_ion_velocity_z = grad_ion_velocity(mach_z, grad_mach_z,
                                                    temperature, grad_temperature)
            ion_vorticity_p =  vc.curl(None,
                                   vector_grad=[grad_ion_velocity_x_p,
                                                grad_ion_velocity_y,
                                                grad_ion_velocity_z])
            ion_vorticity_m =  vc.curl(None,
                                       vector_grad=[grad_ion_velocity_x_m,
                                                    grad_ion_velocity_y,
                                                    grad_ion_velocity_z])

            ion_vorticity_p_smooth = []
            ion_vorticity_m_smooth = []
            for direction in xrange(len(ion_vorticity_p)):
                ion_vorticity_p_smooth.append(gaussian_filter(ion_vorticity_p[direction],
                                                              args.filter_width, mode='reflect'))
                ion_vorticity_m_smooth.append(gaussian_filter(ion_vorticity_m[direction],
                                                              args.filter_width,
                                                              mode='reflect'))
            ion_velocity = ion_velocity_p
            ion_vorticity = ion_vorticity_p
            ion_vorticity_smooth = ion_vorticity_p_smooth

        ## Vector Potential
        ##
        b_field = [fields['B_x'], fields['B_y'], fields['B_z']]
        vector_potential = devore_invert_curl(mesh,
                                              b_field)

        b_field_dynamic = np.array(b_field)
        b_field_dynamic[2] = b_field[2] - args.bias_field_magnitude
        vector_potential_dynamic = devore_invert_curl(mesh,
                                                      b_field_dynamic)
        ## Reference fields
        ##
        (vector_potential_ref,
         b_field_ref,
         b_scalar_potential_ref) = ref_fields(mesh, b_field,
                                              return_scalar_ref=True)
        (vector_potential_dynamic_ref,
         b_field_dynamic_ref,
         b_scalar_potential_dynamic_ref) = ref_fields(mesh, b_field_dynamic,
                                                      return_scalar_ref=True)
        if not just_magnetic:
            (ion_velocity_smooth_ref,
             ion_vorticity_smooth_ref,
             i_scalar_potential_smooth_ref) = ref_fields(mesh, ion_vorticity_smooth,
                                                     return_scalar_ref=True)
            (ion_velocity_ref,
             ion_vorticity_ref,
             i_scalar_potential_ref) = ref_fields(mesh, ion_vorticity,
                                                  return_scalar_ref=True)

            fields = ([b_field[0]] + [b_field[1]] + [b_field[2]] +
                      list(current_smooth) + list(current) +
                      [density] + [temperature] +
                      [density_smooth_plane_normalized] +
                      [temperature_smooth_plane_normalized] +
                      [density] + [temperature] +
                      [density_plane_normalized] +
                      [temperature_plane_normalized] +
                      list(ion_velocity) +
                      list(ion_vorticity_smooth) + list(ion_vorticity) +
                      list(vector_potential) + list(b_field_ref) +
                      list(vector_potential_ref) +
                      list(ion_velocity_ref) +
                      list(ion_vorticity_smooth_ref) +
                      list(ion_vorticity_ref))

            quantity_names = ['B_x', 'B_y', 'B_z',
                              'j_x', 'j_y', 'j_z',
                              'j_raw_x', 'j_raw_y', 'j_raw_z',
                              'n', 'Te',
                              'n_plane_normalized', 'Te_plane_normalized',
                              'n_raw', 'Te_raw',
                              'n_raw_plane_normalized',
                              'Te_raw_plane_normalized',
                              'u_i_x_plus', 'u_i_y', 'u_i_z',
                              'w_i_x_plus', 'w_i_y_plus', 'w_i_z_plus',
                              'w_i_raw_x_plus', 'w_i_raw_y_plus', 'w_i_raw_z_plus',
                              'A_x', 'A_y', 'A_z',
                              'B_ref_x', 'B_ref_y', 'B_ref_z',
                              'A_ref_x', 'A_ref_y', 'A_ref_z',
                              'u_i_ref_x', 'u_i_ref_y', 'u_i_ref_z',
                              'w_i_ref_x', 'w_i_ref_y', 'w_i_ref_z',
                              'w_i_raw_ref_x', 'w_i_raw_ref_y', 'w_i_raw_ref_z']

        else:
            fields = ([b_field[0]] + [b_field[1]] + [b_field[2]] +
                     list(current_smooth) + list(current) +
                     list(b_field_dynamic) +
                     list(b_field_dynamic_ref) +
                     list(vector_potential_dynamic) +
                     list(vector_potential_dynamic_ref) +
                     [b_scalar_potential_ref] +
                     [b_scalar_potential_dynamic_ref])

            quantity_names = ['B_x', 'B_y', 'B_z',
                              'j_x', 'j_y', 'j_z',
                              'j_raw_x', 'j_raw_y', 'j_raw_z',
                              'A_x', 'A_y', 'A_z',
                              'B_ref_x', 'B_ref_y', 'B_ref_z',
                              'A_ref_x', 'A_ref_y', 'A_ref_z',
                              'B_dynamic_x', 'B_dynamic_y', 'B_dynamic_z',
                              'B_dynamic_ref_x', 'B_dynamic_ref_y', 'B_dynamic_ref_z',
                              'A_dynamic_x', 'A_dynamic_y', 'A_dynamic_z',
                              'A_dynamic_ref_x', 'A_dynamic_ref_y', 'A_dynamic_ref_z',
                              'phi_b_ref',
                              'phi_b_dynamic_ref']


        for i, field in enumerate(fields):
            fields[i] = field.astype('float64')

        x, y, z, variables = struc_3d.prepare_for_rectilinear_grid(mesh, fields,
                                                                   quantity_names)
        vtk_file_path = out_dir + args.output_prefix
        struc_3d.write_fields_to_rectilinear_grid(vtk_file_path,
                                                  x, y, z, variables,
                                                  time_point)


def filter_data(data, filter_sigma=None,
                filter_truncate=None):
    r"""
    Filter (with Gaussian).
    """
    if filter_sigma:
        if filter_truncate:
            filtered = ndimage.gaussian_filter(data, filter_sigma,
                                               truncate=filter_truncate)
        else:
            filtered = ndimage.gaussian_filter(data, filter_sigma)
    else:
        filtered = data
    return filtered


def parse_args():
    r"""
    """
    parser = argparse.ArgumentParser(description=("Create VTK files"
                                                  "of canonical quantities"))
    parser.add_argument('--input_path',
                        help='path to input files',
                        default='../output/data_interp_to_rect_grid/')
    parser.add_argument('--input_date',
                        help='time stamp of input files',
                        default='2017-05-10-11-53')
    parser.add_argument('--input_file_text',
                        help='input file name',
                        default='data_interp_to_rect_grid')
    parser.add_argument('--output_prefix',
                        help='prefix of output files',
                        default='canonical_quantities')
    parser.add_argument('--to_cut',
                         help=("number of indicies to cut on"
                               "each side to remove NaNs e.g."
                               "xlow xhigh ylow yhigh zlow zhigh"),
                        nargs=6, type=int,
                        default=[2, -1, 3, -3, 0, -1])
    parser.add_argument('--filter_width',
                        help='width of boxcar filter for derivatives',
                        type=float, default=5)
    parser.add_argument('--interpolate_nan',
                        help='use astropy option to interpolate nans before filtering' +
                        ' may not work well because nans are usually located at edge.',
                        type=bool, default=False)
    parser.add_argument('--density_constant_factor',
                        help='value used for constant density parameters ',
                        type=float,
                        default=1e18)
    parser.add_argument('--bias_field_magnitude',
                        help='magnitude of axial bias magnetic field,' +
                             'used to calculate dynamic field.',
                        type=float,
                        default=0.02)
    parser.add_argument('--time_steps',
                        help='# of time steps', type=int,
                        default=250)
    parser.add_argument('--just_magnetic',
                        help='only calculate magnetic quantities',
                        action='store_true',
                        default=False)
    parser.add_argument('--filter_gaussian',
                        help="run with averaging on all measured quantities,"
                        "used to average measurements for uncertainty analysis",
                        default=False, action='store_true')
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


def normalize_scalar_by_plane(scalar):
    r"""
    Return scalar normalized by the maximum in each z plane.

    Notes
    -----
    To find the max in each plane take nanmax in x and y direction.
    Each nanmax reduces the number of dimensions by 1.
    """
    scalar = np.array(scalar)
    maxes = np.nanmax(np.nanmax(scalar, axis=0), axis=0)
    return scalar / maxes[None, None, :]


def remove_edges_vector(quantity,
                        x_start=2, x_end=None,
                        y_start=0, y_end=-2,
                        z_start=0, z_end=-1):
    r"""
    """
    for index in xrange(len(quantity)):
        quantity[index] = quantity[index][y_start:y_end, x_start:x_end, z_start:z_end]
    return quantity


def remove_edges_scalar(quantity,
                        x_start=2, x_end=None,
                        y_start=0, y_end=-2,
                        z_start=0, z_end=-1):
    r"""
    """
    quantity = quantity[y_start:y_end, x_start:x_end, z_start:z_end]
    return quantity


def boxcar_filter(quantity, width, interpolate_nan=False):
    r"""
    """
    quantity = np.array(quantity)
    nan_indexes = np.where(np.isnan(quantity))
    quantity[nan_indexes] = 0
    boxcar = np.ones((width, width, width))
    boxcar /= boxcar.sum()
    if interpolate_nan:
        return convolve_fft(quantity, boxcar,
                            boundary='wrap',
                            interpolate_nan=True,
                            fftn=scipy.fftpack.fft,
                            ifftn=scipy.fftpack.ifft)
    else:
        return fftconvolve(quantity, boxcar, mode='same')


def grad_ion_velocity(mach, grad_mach, te, grad_te):
    r"""
    """
    factor = q_e/m_i
    term1 = np.sqrt(te*factor)*(grad_mach)
    term2 = mach*factor*grad_te/(2.*np.sqrt(factor*te))
    return term1 + term2




if __name__ == '__main__':
    args = parse_args()
    main(args)
