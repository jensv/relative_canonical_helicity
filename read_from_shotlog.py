# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 14:01:40 2016

@author: Jens von der Linden
"""

import pandas as pd
import sqlite3


def read_and_store_shotlog(shots, database):
    r"""
    """
    for shot in shots:
        mach_settings = read_mach_settings_shotlog(shot)
        store_mach_settings_sql(shot, mach_settings, database)


def read_mach_settings_shotlog(shot):
    r"""
    Return a dictionary with all mach probe settings from the shotlog for a
    specific shot number.
    """
    path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2013-05-02.8.xlsx'
    column = {'mach_insertion': 'Unnamed: 37',
              'mach_x': 'Unnamed: 38',
              'mach_x_read': 'Unnamed: 39',
              'mach_orientation': 'Unnamed: 40',
              'mach_orientation_read': 'Unnamed: 41',
              'mach_y': 'Unnamed: 42',
              'mach_z': 'Unnamed: 43'}
    mach_settings = {}

    shotlog = pd.read_excel(path, sheetname='shot log', index_col=3)

    for key in column.keys():
        mach_settings[key] = shotlog[column[key]][shot]
    return mach_settings


def store_mach_settings_sql(shot, mach_settings, database):
    r"""
    """
    connection = sqlite3.connection()
    cursor = connection.cursor()
    cursor.execute("UPDATE Shots SET mach_x = ?, mach_y = ?, mach_z = ?, " +
                   "mach_orientation = ?  WHERE shot=?;",
                   [mach_settings['mach_x'], mach_settings['mach_y'],
                    mach_settings['mach_z'], shot])
    connection.commit()
    cursor.close()
    connection.close()
