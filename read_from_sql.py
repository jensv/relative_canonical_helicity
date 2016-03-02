# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 09:49:01 2016

@author: Jens von der Linden
"""

import sqlite3

def read_all_rows(condition, database, table_name):
    r"""
    Return all rows from sql table that match condition.
    """
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    cursor = connection.cursor()
    cursor.execute('SELECT * FROM ' + table_name + ' WHERE ' + condition)
    rows = cursor.fetchall()
    cursor.close()
    connection.close()
    return rows


def cursor_with_rows(condition, database, table_name):
    r"""
    Return cursor object which can iterate through rows matching condition.
    """
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    cursor = connection.cursor()
    cursor.execute('SELECT * FROM ' + table_name + ' WHERE ' + condition)
    return cursor, connection


def close(connection, cursor):
    r"""
    Close connection and cursor.
    """
    cursor.close()
    connection.close()
