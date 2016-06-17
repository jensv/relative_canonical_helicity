# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 14:01:40 2016

@author: Jens von der Linden
"""

import pandas as pd
import sqlite3


def read_and_store_shotlog(shots, database, start=False):
    r"""
    """
    path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2013-05-30.1.xlsx'
    shotlog1 = pd.read_excel(path, sheetname='shot log', index_col=3)
    path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2013-05-02.8.xlsx'
    shotlog2 = pd.read_excel(path, sheetname='shot log', index_col=3)
    path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2013-03-08.3.xlsx'
    shotlog3 = pd.read_excel(path, sheetname='shot log', index_col=3)
    path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2012-09-27.xlsx'
    shotlog4 = pd.read_excel(path, sheetname='shot log', index_col=3)
    path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2012-09-14_cumulative.xlsx'
    shotlog5 = pd.read_excel(path, sheetname='shot log', index_col=3)
    shotlogs = [shotlog1, shotlog2, shotlog3, shotlog4, shotlog5]
    for shot in shots:
        print "Writing Shot Settings %i" % shot
        rsx_settings = read_rsx_settings_shotlog(shot, shotlogs)
        store_rsx_settings_sql(shot, rsx_settings, database, start=start)


def read_rsx_settings_shotlog(shot, shotlogs):
    r"""
    Return a dictionary with all mach probe settings from the shotlog for a
    specific shot number.
    """
    rsx_settings = {}
    if 17344 <= shot <= 17622:
        shotlog_num = 0
        path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2013-05-30.1.xlsx'
        column = {'mach_orientation': 'Unnamed: 17',
                  'mach_orientation_read': 'Unnamed: 18',
                  'mach_x': 'Unnamed: 19',
                  'mach_y': 'Unnamed: 20',
                  'mach_y_read': 'Unnamed: 21',
                  'mach_z': 'Unnamed: 22',
                  'bdot10_orientation': 'Unnamed: 23',
                  'bdot10_x': 'Unnamed: 24',
                  'bdot10_y': 'Unnamed: 25',
                  'bdot10_z': 'Unnamed: 26',
                  'tp1_insertion': 'Unnamed: 27',
                  'tp1_x': 'Unnamed: 28',
                  'tp1_y': 'Unnamed: 29',
                  'tp1_z': 'Unnamed: 30',
                  'shotlog': 'Unnamed: 57'}

    elif 16594 <= shot <= 17343:
        shotlog_num = 1
        path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2013-05-02.8.xlsx'
        column = {'bdot3a_orientation': 'Unnamed: 17',
                  'bdot3a_x': 'Unnamed: 18',
                  'bdot3a_y': 'Unnamed: 19',
                  'bdot3a_z': 'Unnamed: 20',
                  'bdot10_orientation': 'Unnamed: 21',
                  'bdot10_x': 'Unnamed: 22',
                  'bdot10_y': 'Unnamed: 23',
                  'bdot10_z': 'Unnamed: 24',
                  'tp1_insertion': 'Unnamed: 25',
                  'tp1_x': 'Unnamed: 26',
                  'tp1_y': 'Unnamed: 27',
                  'tp1_z': 'Unnamed: 28',
                  'mach_insertion': 'Unnamed: 37',
                  'mach_x': 'Unnamed: 38',
                  'mach_x_read': 'Unnamed: 39',
                  'mach_orientation': 'Unnamed: 40',
                  'mach_orientation_read': 'Unnamed: 41',
                  'mach_y': 'Unnamed: 42',
                  'mach_z': 'Unnamed: 43',
                  'shotlog': 'Unnamed: 62'}

    elif 16057 <= shot <= 16578:
        shotlog_num = 2
        path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2013-03-08.3.xlsx'
        column = {'bdot3a_orientation': 'Unnamed: 17',
                  'bdot3a_x': 'Unnamed: 18',
                  'bdot3a_y': 'Unnamed: 19',
                  'bdot3a_z': 'Unnamed: 20',
                  'bdot10_orientation': 'Unnamed: 21',
                  'bdot10_x': 'Unnamed: 22',
                  'bdot10_y': 'Unnamed: 23',
                  'bdot10_z': 'Unnamed: 24',
                  'tp1_insertion': 'Unnamed: 25',
                  'tp1_x': 'Unnamed: 26',
                  'tp1_y': 'Unnamed: 27',
                  'tp1_z': 'Unnamed: 28',
                  'tp2_insertion': 'Unnamed: 37',
                  'tp2_x': 'Unnamed: 38',
                  'tp2_y': 'Unnamed: 39',
                  'tp2_z': 'Unnamed: 40',
                  'shotlog': 'Unnamed: 59'}

    elif 15927 < shot < 16018:
        shotlog_num = 3
        path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2012-09-27.xlsx'
        column = {'bdot3a_orientation': 'Unnamed: 17',
                  'bdot3a_x': 'Unnamed: 18',
                  'bdot3a_y': 'Unnamed: 19',
                  'bdot3a_z': 'Unnamed: 20',
                  'bdot10_orientation': 'Unnamed: 21',
                  'bdot10_x': 'Unnamed: 22',
                  'bdot10_y': 'Unnamed: 23',
                  'bdot10_z': 'Unnamed: 24',
                  'mach_orientation': 'Unnamed: 25',
                  'mach_x': 'Unnamed: 26',
                  'mach_y': 'Unnamed: 27',
                  'mach_z': 'Unnamed: 28',
                  'tp2_insertion': 'Unnamed: 37',
                  'tp2_x': 'Unnamed: 38',
                  'tp2_y': 'Unnamed: 39',
                  'tp2_z': 'Unnamed: 40',
                  'shotlog': 'Unnamed: 59'}

    elif  15249 < shot < 15928:
        shotlog_num = 4
        path = '/Users/vonderlinden2/rsx_drive/RSX/MDS+/Shotlogs/shotlog_2012-09-14_cumulative.xlsx'
        column = {'bdot3a_orientation': 'Unnamed: 17',
                  'bdot3a_x': 'Unnamed: 18',
                  'bdot3a_y': 'Unnamed: 19',
                  'bdot3a_z': 'Unnamed: 20',
                  'bdot10_orientation': 'Unnamed: 21',
                  'bdot10_x': 'Unnamed: 22',
                  'bdot10_y': 'Unnamed: 23',
                  'bdot10_z': 'Unnamed: 24',
                  'tp1_insertion': 'Unnamed: 25',
                  'tp1_x': 'Unnamed: 26',
                  'tp1_y': 'Unnamed: 27',
                  'tp1_z': 'Unnamed: 28',
                  'tp2_insertion': 'Unnamed: 37',
                  'tp2_x': 'Unnamed: 38',
                  'tp2_y': 'Unnamed: 39',
                  'tp2_z': 'Unnamed: 40',
                  'shotlog': 'Unnamed: 59'}
    else:
        print 'No column dictionary defined for %i' % shot
        return None

    if shot in shotlogs[shotlog_num]['Unnamed: 0'].index:
        for key in column.keys():
            rsx_settings[key] = shotlogs[shotlog_num][column[key]][shot]
        rsx_settings['exists_in_shotlog'] = True
    else:
        print 'Shot %i not in excel spreadsheet.' % shot
        rsx_settings['exists_in_shotlog'] = False
    return rsx_settings


def store_rsx_settings_sql(shot, rsx_settings, database, start=False):
    r"""
    Write contents of one row (shot) from shotlog sheet to sql database.
    """
    connection = sqlite3.connect(database)
    cursor = connection.cursor()
    if start:
        if rsx_settings:
            key_string = '('
            keys = rsx_settings.keys()
            for key in keys[:-1]:
                key_string += ":" + key + ", "
            key_string += ":" + keys[-1] + ", :shot)"
            rsx_settings.update({'shot': shot})
            cursor.execute("INSERT INTO Shots ( " + str(keys + ['shot'])[1:-1].replace("'", "") +
                           ") VALUES " + key_string + ";",
                           rsx_settings)
        else:
            cursor.execute("INSERT INTO Shots (shot, exists_in_shotlog) " +
                           "VALUES (:shot, :exists);", {'shot': shot, 'exists': False})
    else:
        if rsx_settings:
            for key in rsx_settings.keys():
                cursor.execute("UPDATE Shots SET " + key + " = :value " +
                               "WHERE shot = :shot;",
                               {'value': rsx_settings[key], 'shot': shot})
        else:
            cursor.execute("UPDATE Shots SET exists_in_shotlog = :exists " +
                           "WHERE shot = :shot;", {'shot': shot, 'exists': False})
    connection.commit()
    cursor.close()
    connection.close()
