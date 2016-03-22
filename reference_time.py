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


def determine_reference_times_all_shots(shots, database, plot=True,
                                        time_range=[1.5, 2.5], start=False):
    r"""
    Determine relative times for each shot using Jason's script, store in sql
    database.
    """
    for shot in shots:
        connection = sqlite3.connect(database)
        cursor = connection.cursor()
        cursor.execute("SELECT fiducial_a_node, fiducial_b_node, " +
                       "bias_current_node FROM Shots WHERE shot = :shot;",
                       {'shot': shot})
        fiducial_a_node, fiducial_b_node, bias_current_node = cursor.fetchone()
        cursor.close()
        connection.close()

        print 'shot:', shot
        determine_time, existence = shot_exists(shot, fiducial_a_node,
                                                fiducial_b_node,
                                                bias_current_node)
        if determine_time:
            times, success = determine_times_from_idl(shot)
            store_times_in_sql(shot, times, database, start=start)
            if plot:
                raw_signals = retrieve_fiducial_signals(shot, fiducial_a_node,
                                                        fiducial_b_node,
                                                        bias_current_node)
                save_as = '../output/' + str(shot) + '.png'
                if success:
                    pass
                    #plot_signals(raw_signals,
                    #             zero_times=[times['phase_zero'], times['ramp']],
                    #             legend=True, show=False, time_range=time_range,
                    #             save_as=save_as)
                else:
                    pass
                    #plot_signals(raw_signals,
                    #             legend=True, show=False, time_range=time_range,
                    #             save_as=save_as)
        store_existence_in_sql(shot, database, existence)


def shot_exists(shot, fiducial_a_node, fiducial_b_node, bias_current_node):
    r"""
    Tries to read data from fiducial and current monitor nodes of rsx MDSplus tree.
    If opening of tree or reading from nodes errors returns False otherwise True.
    """
    exists_in_mdsplus = True
    fiducial_a_signal_exists = True
    fiducial_b_signal_exists = True
    bias_current_signal_exists = True
    try:
        rsx_tree = mds.Tree('rsx', shot)
    except:
        print '%i does not exist' % shot
        exists_in_mdsplus = False
    if exists_in_mdsplus:
        try:
            node = rsx_tree.getNode(fiducial_a_node).getData()
            assert not len(node.getValue().getShape()) == 0
            assert node.getValue().size >= 10000
        except:
            fiducial_a_signal_exists = False
        try:
            node = rsx_tree.getNode(fiducial_b_node).getData()
            assert not len(node.getValue().getShape()) == 0
            assert node.getValue().size >= 10000
        except:
            fiducial_b_signal_exists = False
        try:
            node = rsx_tree.getNode(bias_current_node).getData()
            assert not len(node.getValue().getShape()) == 0
            assert node.getValue().size >= 10000
        except:
            bias_current_signal_exists = False
    determine_time = True  if (exists_in_mdsplus and
                               fiducial_a_signal_exists and
                               bias_current_signal_exists) else False
    existence = {'exists_in_mdsplus': exists_in_mdsplus,
                 'fiducial_a_signal_exists': fiducial_a_signal_exists,
                 'fiducial_b_signal_exists': fiducial_b_signal_exists,
                 'bias_current_signal_exists': bias_current_signal_exists}
    return (determine_time, existence)


def determine_times_from_idl(shot,
                             idl_path='/Applications/exelis/IDL85/bin/idl'):
    r"""
    Run relative time script and return dictionary of ouputs.
    """
    success = True
    times = {}
    try:
        idl = pidly.IDL(idl_path)
        idl.pro('pro00100, '+str(shot)+', trig_index, trig_time, period;')
        phase_zero_time = float(idl.trig_time)*1e-3
        phase_zero_index = int(idl.trig_index)
        period = float(idl.period)*1e-3
        idl.close()
        times.update({'phase_zero': phase_zero_time,
                      'phase_zero_index': phase_zero_index,
                      'period': period,
                      'phase_reference_time_idl_code_succeeded': True})
    except:
        times.update({'phase_zero': None,
                      'phase_zero_index': None,
                      'period': None,
                      'phase_reference_time_idl_code_succeeded': False})
        success = False
    try:
        idl = pidly.IDL(idl_path)
        idl.pro('pro00100, '+str(shot)+', trig_index, trig_time, period, current_rise=1;')
        ramp_time = float(idl.trig_time)*1e-3
        ramp_index = int(idl.trig_index)
        idl.close()
        times.update({'ramp': ramp_time,
                      'ramp_index': ramp_index,
                      'ramp_reference_time_idl_code_succeeded': True})
    except:
        times.update({'ramp': None,
                      'ramp_index': None,
                      'ramp_reference_time_idl_code_succeeded': False})
        success = False
    return times, success


def store_times_in_sql(shot, times, database, start=False):
    r"""
    Store times in sql database.
    """
    connection = sqlite3.connect(database)
    cursor = connection.cursor()
    if start:
        insert_statement = ("INSERT INTO Shots (shot, existence, " +
                            "zero_phase_time, zero_phase_index, period, "
                            "ramp_time, ramp_index) " +
                            "VALUES (?, ?, ?, ?, ?, ?, ?);")
        cursor.execute(insert_statement, [shot, True, times['phase_zero'],
                       times['phase_zero_index'],
                       times['period'], times['ramp'],
                       times['ramp_index']])
    else:
        update_statement = ("UPDATE Shots SET zero_phase_time = :phase_zero, " +
                            "zero_phase_index = :phase_zero_index, " +
                            "period = :period, ramp_time = :ramp_time, " +
                            "ramp_index = :ramp_index, " +
                            "phase_reference_time_idl_code_succeeded = :phase_reference_time_idl_code_succeeded, " +
                            "ramp_reference_time_idl_code_succeeded = :ramp_reference_time_idl_code_succeeded, " +
                            "WHERE shot = :shot;")
        cursor.execute(update_statement, {'shot': shot,
                                          'phase_reference_time_idl_code_succeeded': times['phase_reference_time_idl_code_succeeded'],
                                          'ramp_reference_time_idl_code_succeeded' : times['ramp_reference_time_idl_code_succeeded'],
                                          'phase_zero': times['phase_zero'],
                                          'phase_zero_index':
                                              times['phase_zero_index'],
                                          'period': times['period'],
                                          'ramp_time': times['ramp'],
                                          'ramp_index': times['ramp_index']})
    connection.commit()
    cursor.close()
    connection.close()


def store_existence_in_sql(shot, database, existence):
    r"""
    Store shot with existense set to False.
    """
    connection = sqlite3.connect(database)
    cursor = connection.cursor()
    insert_statement = ("UPDATE Shots SET exists_in_mdsplus = :exists_in_mdsplus, " +
                        "fiducial_a_signal_exists = :fiducial_a_signal_exists, " +
                        "fiducial_b_signal_exists = :fiducial_b_signal_exists, " +
                        "bias_current_signal_exists = :bias_current_signal_exists " +
                        "WHERE shot = :shot;")
    existence.update({'shot': shot})
    cursor.execute(insert_statement, existence)
    connection.commit()
    cursor.close()
    connection.close()


def retrieve_fiducial_signals(shot, fiducial_a_node, fiducial_b_node,
                              bias_current_node):
    r"""
    Read MDS+ tree nodes storing fiducial probes and return signals as dictionary.
    """
    rsx_tree = mds.Tree('rsx', shot)
    fiducial_a_node = rsx_tree.getNode(fiducial_a_node)
    fiducial_b_node = rsx_tree.getNode(fiducial_b_node)
    bias_current_node = rsx_tree.getNode(bias_current_node)
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
            plt.axvline(zero_time*1e3)
    if show:
        plt.show()
    if save_as:
        plt.savefig(save_as)
    plt.close()
