# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:52:04 2016

@author: Jens von der Linden
"""

import sqlite3

connection = sqlite3.connect('shots.db')

cursor = connection.cursor()

cursor.execute("CREATE TABLE Shots(shot INTEGER, exists_in_shotlog BOOLEAN, " +
               "exists_in_mdsplus BOOLEAN, " +
               "campaign STRING, fiducial_a_signal_exists BOOLEAN, " +
               "fiducial_b_signal_exists BOOLEAN, bias_current_signal_exists, " +
               "bias_current_on BOOLEAN, "+
               "bias_current_pre_crowbar_regular BOOLEAN, " +
               "fiducial_ocillates BOOLEAN, " +
               "zero_phase_time REAL, zero_phase_index INTEGER, " +
               "period FLOAT, ramp_time FLOAT, ramp_index INTEGER, " +
               "mach_insertion FLOAT, " +
               "mach_signals_exist BOOLEAN, mach_oscillates BOOLEAN, " +
               "mach_z_read FLOAT, mach_x FLOAT, mach_y_read FLOAT, mach_y FLOAT, mach_z FLOAT, " +
               "mach_orientation FLOAT, mach_orientation_read FLOAT, mach_r_node STRING, " +
               "mach_l_node STRING, mach_l_node_2 STRING, " +
               "mach_r_termination_ohm FLOAT, mach_l_termination_ohm FLOAT, " +
               "mach_r_monitor_pearson_model INTEGER, mach_l_monitor_pearson_model INTEGER, " +
               "mach_l_2_termination_ohm FLOAT, mach_l_2_monitor_pearson_model INTEGER, " +
               "mach_r_loops INTEGER, mach_l_loops INTEGER, mach_l_2_loops INTEGER, " +
               "mach_r_monitor_volts_per_amp FLOAT, mach_l_monitor_volts_per_amp FLOAT," +
               "mach_l_2_monitor_volts_per_amp FLOAT, " +
               "tp1_signals_exist BOOLEAN, " +
               "tp2_signals_exist BOOLEAN, " +
               "tp1_insertion FLOAT, tp1_x Float, tp1_y FLOAT, tp1_z FLOAT, " +
               "tp2_insertion FLOAT, tp2_x FLOAT, tp2_y FLOAT, tp2_z FLOAT, " +
               "bdot3a_signals_exist BOOLEAN, " +
               "bdot3a_x FLOAT, bdot3a_y FLOAT, bdot3a_z FLOAT, bdot3a_orientation FLOAT, " +
               "bdot10_signals_exist BOOLEAN, " +
               "bdot10_x FLOAT, bdot10_y FLOAT, bdot10_z FLOAT, " +
               "bdot10_orientation FLOAT, shotlog STRING);")
connection.commit()

cursor.close()
connection.close()
