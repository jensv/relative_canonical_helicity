# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 17:45:11 2016

@author: Jens von der Linden
"""

import numpy as np
import matplotlib.pyplot as plt
import pidly
import MDSplus as mds
import sqlite3

def determine_relative_times_all_shots(shots, plot=False):
    r"""
    Determine relative times for each shot using Jason's script, store in sql database.
    """
    for shot in shots:
        if shot_exists(shot):
            times = determine_times_from_idl(shot)
            store_in_sql(shot, times)
            if plot:
                raw_signals = retrieve_signals(shot)
                plot(raw_signals, time_range, log_sacle, times)
        else:
            store_nonexistence_in_sql(shot)


def shot_exists(shot):
    r"""
    Tries to open MDSplus tree, if opening errors returns False otherwise True.
    """
    try:
        mds.Tree('rsx', shot)
    except:
        return False
    return True


def determine_times_from_idl(shot,
                             idl_path='/Applications/exelis/IDL85/bin/idl'):
    r"""
    Run relative time script and return dictionary of ouputs.
    """
    idl = pidly.IDL(idl_path)
    idl.pro('pro00100, 16670, trig_index, trig_time, period;')
    phase_zero_time = idl.trig_time[0]
    phase_zero_index = idl.trig_index[0]
    period = idl.period[0]
    idl.pro('pro00100, 16670, trig_index, trig_time, period, current_rise=1;')
    ramp_time = idl.trig_time[0]
    ramp_index = idl.trig_index[0]
    idl.close()
    times = {'phase_zero': phase_zero_time,
             'phase_zero_index': phase_zero_index,
             'ramp': ramp_time,
             'ramp_index': ramp_index}
    return times


def store_times_in_sql():
    r"""
    """
    pass


def retrieve_fiducial_signals(shot):
    r"""
    Read MDS+ tree nodes storing fiducial probes and return signals as dictionary.
    """
    rsx_tree = mds.Tree('rsx', shot)
    fiducial_a_node = rsx_tree.getNode('\j_002_000')
    fiducial_b_node = rsx_tree.getNode('\j_002_001')
    bias_current_node = rsx_tree.getNode('\j_002_004')
    fiducial_a_data = fiducial_a_node.getData()
    fiducial_b_data = fiducial_b_node.getData()
    bias_current_data = bias_current_node.getData()
    fiducial_a_raw = np.asarray(fiducial_a_data.getValue())
    fiducial_b_raw = np.asarray(fiducial_b_data.getValue())
    bias_current_raw = np.asarray(bias_current_data.getValue())

    fiducial_a_time = np.asarray(fiducial_a_data.getDimensions()[0])*1e-3
    fiducial_b_time = np.asarray(fiducial_b_data.getDimensions()[0])*1e-3
    bias_current_time = np.asarray(bias_current_data.getDimensions()[0])*1e-3
    assert np.allclose(fiducial_a_time,
                       fiducial_b_time) and np.allclose(fiducial_a_time,
                                                        bias_current_time)
    raw_signals = {'fid_a': fiducial_a_raw,
                   'fid_b': fiducial_b_raw,
                   'bias': bias_current_raw,
                   'time': fiducial_a_time}
    return raw_signals


def plot_signals(raw_signals, log_scale=True, time_range=None,
                 fid_b=False, zero_times=None, legend=False, show=True,
                 save_as=None):
    r"""
    """
    plt.plot(raw_signals['time']*1e3, raw_signals['fid_a'],
             label='fiducial probe a')
    if fid_b:
        plt.plot(raw_signals['time']*1e3, raw_signals['fid_b'],
                 label='fiducial probe b')
    plt.plot(raw_signals['time']*1e3, raw_signals['bias'],
             label='bias current')
    plt.ylabel('raw Voltage [V]')
    plt.xlabel('time [ms]')
    if log_scale:
        plt.yscale('log')
    if time_range:
        plt.xlim(time_range)
    if legend:
        plt.legend()
    if zero_times:
        np.asarray(zero_times)
        for zero_time in zero_times:
            plt.axvline(zero_time)
    if show:
        plt.show()
    if save_as:
        plt.savefig(save_as)