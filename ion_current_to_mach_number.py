# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 18:20:47 2016

@author: Jens von der Linden
"""

import numpy as np
import MDSplus as mds


def mach_number(shot, droop_factor_l=0.8e3, droop_factor_r=0.7e3,
                mach_calibration_factor=0.45, zero_signal_points=2999,
                tree_name='rsx', mach_r_node_name='\j_008_015',
                mach_l_node_name='\j_008_014'):
    r"""
    Top level mach probe analysis function.
    """
    (r_raw, r_time,
     l_raw, l_time) = read_isat(shot, tree_name, mach_r_node_name,
                                mach_l_node_name)
    assert np.allclose(r_time, l_time), 'Sample times for left and right are not identical.'
    r_raw_no_offset, r_background_std = remove_offset(r_raw,
                                                      zero_signal_points)
    l_raw_no_offset, l_background_std = remove_offset(l_raw,
                                                      zero_signal_points)
    r_isat = correct_droop(r_raw_no_offset, r_time, droop_factor_r)
    l_isat = correct_droop(l_raw_no_offset, l_time, droop_factor_l)
    mach = determine_mach(r_isat, l_isat, mach_calibration_factor)
    return mach, r_time, r_background_std, l_background_std


def read_isat(shot, tree_name, mach_r_node_name, mach_l_node_name):
    r"""
    Read saturation current measurements from MDS+ database.
    """
    rsx_tree = mds.Tree(tree_name, shot)
    mach_r_node = rsx_tree.getNode(mach_r_node_name)
    mach_l_node = rsx_tree.getNode(mach_l_node_name)
    mach_r_data = mach_r_node.getData()
    mach_l_data = mach_l_node.getData()
    mach_r_raw = np.asarray(mach_r_data.getValue())
    mach_l_raw = np.asarray(mach_l_data.getValue())
    mach_r_time = np.asarray(mach_r_data.getDimensions()[0])*1e-3
    mach_l_time = np.asarray(mach_l_data.getDimensions()[0])*1e-3
    return mach_r_raw, mach_r_time, mach_l_raw, mach_l_time


def remove_offset(raw_signal, zero_signal_points):
    r"""
    Shift the signal so that it is centered around zero.
    """
    offset = raw_signal[:zero_signal_points].mean()
    raw_signal_no_offset = raw_signal - offset
    background_std = raw_signal[:zero_signal_points].std()
    msg = 'zero_signal_points includes measurements where signal exceeds three times the std.'
    assert np.sum(raw_signal_no_offset > 3. * background_std), msg
    return raw_signal_no_offset, background_std


def correct_droop(uncorrected_signal, signal_time, calibration_factor):
    r"""
    Apply droop correction.
    """
    dt = signal_time[1] - signal_time[0]
    msg = 'dt is not uniform along signal time'
    assert np.allclose(np.diff(signal_time), dt), msg
    signal_integral = np.cumsum(uncorrected_signal)*dt
    correction_term = signal_integral * calibration_factor
    droop_corrected_signal = uncorrected_signal + correction_term
    return droop_corrected_signal


def determine_mach(i_sat_r, i_sat_l, calibration_factor):
    r"""
    Take log of ratio of saturation currents and apply Mach calibration factor.
    """
    mach = np.log(i_sat_l / i_sat_r)*calibration_factor
    return mach
