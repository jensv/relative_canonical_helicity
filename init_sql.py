# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:52:04 2016

@author: Jens von der Linden
"""

import sqlite3

connection = sqlite3.connect('shots.db')

cursor = connection.cursor()

cursor.execute("CREATE TABLE Shots(shot INTEGER, existence BOOLEAN, " +
               "campaign STRING, bias_current_on BOOLEAN, "+
               "bias_current_pre_crowbar_regular BOOLEAN, "+
               "fiducial_ocillates BOOLEAN, " +
               "zero_phase_time REAL, zero_phase_index INTEGER, " +
               "period FLOAT, ramp_time FLOAT, ramp_index INTEGER, " +
               "mach_aquired BOOLEAN, mach_oscillates BOOLEAN, "
               "mach_x FLOAT, mach_y FLOAT, mach_z FLOAT, " +
               "mach_orientation FLOAT);")
connection.commit()

cursor.close()
connection.close()

