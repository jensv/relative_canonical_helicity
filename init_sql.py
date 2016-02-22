# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:52:04 2016

@author: Jens von der Linden
"""

import sqlite3

connection = sqlite3.connect('../output/relative_times.db')

cursor = connection.cursor()

cursor.execute("CREATE TABLE RelativeTimes(shot INTEGER, existence BOOLEAN, " +
               "zero_phase_time REAL, zero_phase_index INTEGER, " +
               "period FLOAT, ramp_time FLOAT, ramp_index INTEGER);")
connection.commit()

cursor.close()
connection.close()

