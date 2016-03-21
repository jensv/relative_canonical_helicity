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
from read_from_shotlog import read_and_store_shotlog
from reference_time import determine_reference_times_all_shots


def build_database(database='shots.db', start=15253, end=17623):
    r"""
    Call functions that populate database.
    """
    read_and_store_shotlog(range(start, end), database=database, start=True)
    print "Classifying by campaign"
    add_campaigns(database)
    add_manual_entries(database)
    add_reference_times(database)


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
    if 15928 <= shot <= 16017:
        campaigns.append('mach_line_campaign_1')
    elif 16596 <= shot <= 16627:
        campaigns.append('mach_probe_point_campaign_1')
    elif 16628 <= shot <= 16646:
        campaigns.append('mach_probe_point_campaign_2')
    elif 16647 <= shot <= 16669:
        campaigns.append('mach_probe_point_campaign_3')
    elif 16670 <= shot <= 16689:
        campaigns.append('mach_probe_point_campaign_4')
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
    determine_reference_times_all_shots(shots, database, start=False, plot=False)


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
            manual_settings['mach_r_node'] = 'j_008_011'
            manual_settings['mach_l_node'] = 'j_008_010'
            manual_settings['mach_l_2_node'] = 'j_002_003'
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
            manual_settings['mach_l_2_monitor_volts_per_amp'] = 1
        if 16579 <= shot:
            manual_settings['mach_r_node'] = 'j_008_015'
            manual_settings['mach_l_node'] = 'j_008_014'
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
