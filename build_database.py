# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 08:48:53 2016

@author: Jens von der Linden
"""

import sys

sys.path.append('../../time_alignment/source')
sys.path.append('../../read_from_shotlog/')

import sqlite3
import numpy as np
import MDSplus as mds
import scipy.signal as signal
import scipy.fftpack as fftpack
import scipy.integrate as integrate
from read_from_shotlog import read_and_store_shotlog
from reference_time import determine_reference_times_all_shots


def build_database(database='shots.db', start=15253, end=17623):
    r"""
    Call functions that populate database.
    """
    heuristics_params = {'pre_ramp_slice': slice(20000, 30000),
                         'ramp_slice': slice(30000, None),
                         'crowbar_slice': slice(5000, None),
                         'peak_smooth_window': 10,
                         'diff_smooth_window': 750,
                         'offset': 5e-6,
                         'pre_crowbar_sd_slice': slice(35500, 60000),
                         'sd_slice': slice(35500, 42300),
                         'mean_gyration_freq': 58e3,
                         'std_gyration_freq': 8.4e3,
                         'ddof': 1}

    read_and_store_shotlog(range(start, end), database=database, start=True)
    print "Classifying by campaign"
    add_campaigns(database)
    add_manual_entries(database)
    add_reference_times(database)
    add_shot_quality_heuristics(database, heuristics_params)
    add_mach_quality_heuristics(database)


def add_mach_quality_heuristics(database):
    r"""
    Tries to read data from mach nodes of rsx MDSplus tree.
    If opening of tree or reading from nodes errors returns False otherwise True.
    """
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM Shots WHERE exists_in_mdsplus=1 AND "
                   "(campaigns = 'mach_probe_plane_campaign_1' OR "
                   "campaigns = 'mach_probe_plane_campaign_2' OR "
                   "campaigns = 'mach_probe_line_campaign_1' OR "
                   "campaigns = 'mach_probe_point_campaign_1' OR "
                   "campaigns = 'mach_probe_point_campaign_2' OR "
                   "campaigns = 'mach_probe_point_campaign_3' OR "
                   "campaigns = 'mach_probe_point_campaign_4' OR "
                   "campaigns = 'mach_probe_point_campaign_5')")
    rows = cursor.fetchall()

    for row in rows:
        mach_signals_exist = True
        try:
            rsx_tree = mds.Tree('rsx', row['shot'])
        except:
            print 'MDS+ error: %i does not exist' % row['shot']

        try:
            l_node = rsx_tree.getNode(row['mach_l_node']).getData()
            assert not len(l_node.getValue().getShape()) == 0, str(row['shot']) + ' l_node_data length zero.'
            assert l_node.getValue().size >= 10000, str(row['shot']) + ' l_node_data short.'

            r_node = rsx_tree.getNode(row['mach_r_node']).getData()
            assert not len(r_node.getValue().getShape()) == 0, str(row['shot']) + ' r_node_data length zero.'
            assert r_node.getValue().size >= 10000, str(row['shot']) + ' r_node_data short.'
        except:
            mach_signals_exist = False

        print 'updating' + str(row['shot']) + ' : ' +  str(mach_signals_exist)
        cursor.execute('UPDATE Shots SET mach_signals_exist = :signals_exist '
                       'WHERE shot = :shot',
                       {'signals_exist': mach_signals_exist,
                        'shot': row['shot']})
    connection.commit()
    cursor.close()
    connection.close()


def add_campaigns(database):
    r"""
    """
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM Shots")
    rows = cursor.fetchall()
    for row in rows:
        campaigns = assign_campaigns(row['shot'])
        cursor.execute('UPDATE Shots SET campaigns = ? WHERE shot = ?',
                       [campaigns, row['shot']])
    connection.commit()
    cursor.close()
    connection.close()


def assign_campaigns(shot):
    r"""
    """
    campaigns = []
    if 15928 <= shot <= 15941:
        campaigns.append('mach_probe_point_campaign_1')
    if 15942 <= shot <= 16017:
        campaigns.append('mach_probe_line_campaign_1')
    elif 16596 <= shot <= 16627:
        campaigns.append('mach_probe_point_campaign_2')
    elif 16628 <= shot <= 16646:
        campaigns.append('mach_probe_point_campaign_3')
    elif 16647 <= shot <= 16669:
        campaigns.append('mach_probe_point_campaign_4')
    elif 16670 <= shot <= 16689:
        campaigns.append('mach_probe_point_campaign_5')
    elif 16690 <= shot <= 17343:
        campaigns.append('mach_probe_plane_campaign_1')
    elif 17344 <= shot <= 17622:
        campaigns.append('mach_probe_plane_campaign_2')
    if 15253 <= shot <= 15927:
        campaigns.append('bdot_tp_plane_campaign_1')
    elif 16057 <= shot <= 16578:
        campaigns.append('bdot_tp_plane_campaign_2')
    campaigns_string = ''
    for campaign in campaigns:
        campaigns_string += campaign + ', '
    campaigns_string = campaigns_string[:-2]
    return campaigns_string


def add_reference_times(database):
    r"""
    """
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM Shots")
    rows = cursor.fetchall()
    shots = [row['shot'] for row in rows]
    determine_reference_times_all_shots(shots, database, start=False,
                                        plot=False)


def add_manual_entries(database):
    r"""
    """
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM Shots")
    rows = cursor.fetchall()
    for row in rows:
        shot = row['shot']
        print "Add manual entries to shot %i" % shot
        manual_settings = {'fiducial_a_node': '\j_002_000',
                           'fiducial_b_node': '\j_002_001',
                           'bias_current_node': '\j_002_004'}
        if 15928 <= shot <= 16017:
            manual_settings['mach_r_node'] = '\j_008_011'
            manual_settings['mach_l_node'] = '\j_008_010'
            manual_settings['mach_l_2_node'] = '\j_002_003'
            manual_settings['mach_r_termination_ohm'] = 95e3
            manual_settings['mach_l_termination_ohm'] = 95e3
            manual_settings['mach_r_monitor_pearson_model'] = '4100'
            manual_settings['mach_l_monitor_pearson_model'] = '4100'
            manual_settings['mach_l_2_termination_ohm'] = 95e3
            manual_settings['mach_l_2_monitor_pearson_model'] = '4110'
            manual_settings['mach_r_loops'] = 1
            manual_settings['mach_l_loops'] = 1
            manual_settings['mach_l_2_loops'] = 10
            manual_settings['mach_r_monitor_volts_per_amp'] = 1
            manual_settings['mach_l_monitor_volts_per_amp'] = 1
            manual_settings['mach_l_2_monitor_volts_per_amp'] = 0.1
        if shot == 15948:
            manual_settings['mach_l_loops'] = 10
        if 16579 <= shot:
            manual_settings['mach_r_node'] = '\j_008_015'
            manual_settings['mach_l_node'] = '\j_008_014'
            manual_settings['mach_r_termination_ohm'] = 95e3
            manual_settings['mach_l_termination_ohm'] = 95e3
            manual_settings['mach_r_monitor_pearson_model'] = '4100'
            manual_settings['mach_l_monitor_pearson_model'] = '4100'
            manual_settings['mach_r_loops'] = 1
            manual_settings['mach_l_loops'] = 1
            manual_settings['mach_r_monitor_volts_per_amp'] = 1
            manual_settings['mach_l_monitor_volts_per_amp'] = 1
        for key in manual_settings.keys():
            cursor.execute("UPDATE Shots SET " + key + " = :value " +
                           "WHERE shot = :shot;",
                           {'value': manual_settings[key], 'shot': shot})
    connection.commit()
    cursor.close()
    connection.close()


def add_shot_quality_heuristics(database, heuristics_params):
    r"""
    """
    add_bias_current_heuristics(database, **heuristics_params)
    add_fiducial_signal_heuristics(database, **heuristics_params)


def add_bias_current_heuristics(database, pre_ramp_slice, ramp_slice,
                                crowbar_slice, peak_smooth_window,
                                diff_smooth_window, offset, ddof, **kwargs):
    r"""
    """
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM Shots WHERE exists_in_shotlog = 1 AND " +
                   "exists_in_mdsplus = 1 AND fiducial_a_signal_exists = 1 AND " +
                   "bias_current_signal_exists = 1;")
    rows = cursor.fetchall()
    for row in rows:
        print "Add bias current heuristics to shot %i" % row['shot']
        bias_current_node_name = row['bias_current_node']
        tree = mds.Tree('rsx', row['shot'])
        bias_current_node = tree.getNode(bias_current_node_name)
        bias_current_data = bias_current_node.getData()
        bias_current = np.asarray(bias_current_data.getValue())
        bias_current_time =  np.asarray(bias_current_data.getDimensions()[0])*1e-3
        std = get_bias_current_std(bias_current_time, bias_current,
                                   pre_ramp_slice, ddof)
        peak_current = get_bias_current_peak(bias_current_time, bias_current,
                                             ramp_slice, peak_smooth_window)
        crowbar_time = get_crowbar_time(bias_current_time, bias_current,
                                        ramp_slice,
                                        crowbar_slice, offset,
                                        diff_smooth_window)
        std = std.tolist()
        peak_current = peak_current.tolist()
        if crowbar_time:
            crowbar_time = crowbar_time.tolist()
        cursor.execute("UPDATE Shots SET bias_current_pre_ramp_std = :std, " +
                       "bias_current_peak = :peak_current, " +
                       "bias_current_crowbar_time = :crowbar_time " +
                       "WHERE shot=:shot",
                       {'shot': row['shot'], 'std': std, 'peak_current':
                        peak_current, 'crowbar_time': crowbar_time})
    connection.commit()
    cursor.close()
    connection.close()


def get_bias_current_std(bias_current_time, bias_current, pre_ramp_slice,
                         ddof):
    r"""
    """
    return bias_current[pre_ramp_slice].std(ddof=ddof)


def get_bias_current_peak(bias_current_time, bias_current, ramp_slice,
                          window_length):
    r"""
    """
    window = np.ones(window_length) / window_length
    smoothed = np.convolve(window, bias_current[ramp_slice], mode='valid')
    return smoothed.max()


def get_crowbar_time(bias_current_time, bias_current, ramp_slice,
                     crowbar_slice, offset, window_length):
    r"""
    """
    window = np.ones(window_length) / window_length
    smoothed = np.convolve(window, bias_current[ramp_slice], mode='same')
    derivative = np.gradient(smoothed)[crowbar_slice] + offset
    try:
        zero_crossing = np.where(np.diff(np.sign(derivative)))[0][0]
        crowbar_time = bias_current_time[ramp_slice][crowbar_slice][zero_crossing]
    except:
        crowbar_time = None
    return crowbar_time


def add_fiducial_signal_heuristics(database, pre_crowbar_sd_slice, sd_slice,
                                   std_gyration_freq, mean_gyration_freq,
                                   **kwargs):
    r"""
    """
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM Shots WHERE exists_in_shotlog = 1 AND " +
                   "exists_in_mdsplus = 1 AND fiducial_a_signal_exists = 1 AND " +
                   "bias_current_signal_exists = 1 AND period IS NOT NULL;")
    rows = cursor.fetchall()
    for row in rows:
        print "Add fiducial signal heuristics to shot %i" % row['shot']
        period = row['period']
        zero_phase_index = row['zero_phase_index']
        gyration_freq = 1. / row['period']
        tree = mds.Tree('rsx', row['shot'])
        fiducial_node_name = row['fiducial_a_node']
        fiducial_node = tree.getNode(fiducial_node_name)
        fiducial_data = fiducial_node.getData()
        fiducial = np.asarray(fiducial_data.getValue())
        fiducial_time =  np.asarray(fiducial_data.getDimensions()[0])*1e-3
        time_step = fiducial_time[1] - fiducial_time[0]
        pre_crowbar_gyration_sd = get_gyration_sd(fiducial_time, fiducial,
                                                  sd_slice, gyration_freq,
                                                  std_gyration_freq)
        gyration_sd = get_gyration_sd(fiducial_time, fiducial,
                                      pre_crowbar_sd_slice, gyration_freq,
                                      std_gyration_freq)
        n_std = np.ceil((gyration_freq - mean_gyration_freq) / std_gyration_freq)

        ### determine uncalibrated gyration amplitude
        start_index = zero_phase_index
        end_index = zero_phase_index + np.ceil(period/time_step)*1.25
        zero = np.mean(fiducial[:500])
        integral = cumtrapz(fiducial[start_index:end_index]-zero, dx=time_step, initial=0)
        gyration_amplitude = np.min(integral)/2.

        cursor.execute("UPDATE Shots SET " +
                       "fiducial_pre_crowbar_gyration_spectral_density = :pre_sd, " +
                       "fiducial_gyration_spectral_density = :sd, " +
                       "period_within_n_std = :n_std, " +
                       "uncalibrated_integrated_fiducial_a_gyration_amplitude = :gyration_amplitude "
                       "WHERE shot = :shot;",
                       {'shot': row['shot'], 'pre_sd': pre_crowbar_gyration_sd.tolist(),
                        'sd': gyration_sd.tolist(), 'n_std': n_std.tolist(),
                        'gyration_amplitude': gyration_amplitude})
    connection.commit()
    cursor.close()
    connection.close()


def get_gyration_sd(fiducial_time, fiducial,
                    sd_slice, gyration_freq, std_gyration_freq):
    r"""
    """
    window = fiducial[sd_slice]
    fs = 1/(fiducial_time[1] - fiducial_time[0])
    freqs, periodogram = signal.periodogram(window, fs)
    range_of_interest = np.where(np.logical_and(gyration_freq - std_gyration_freq <= freqs,
                                                freqs <= gyration_freq + std_gyration_freq))

    sd = integrate.cumtrapz(periodogram[range_of_interest], freqs[range_of_interest])[0]
    return sd
