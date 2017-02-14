# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 10:43:31 2016

@author: Jens von der Linden
"""


def maximum_b_extent():
    extent_0249_dict = {'x_min': -0.028, 'x_max': 0.025,
                        'y_min': -0.043, 'y_max': 0.039,
                        'z_min': 0.249, 'z_max': 0.249}
    extent_0302_dict = {'x_min': -0.022, 'x_max': 0.021,
                        'y_min': -0.038, 'y_max': 0.04,
                        'z_min': 0.302, 'z_max': 0.302}
    extent_0357_dict = {'x_min': -0.041, 'x_max': 0.030,
                        'y_min': -0.019, 'y_max': 0.0255,
                        'z_min': 0.357, 'z_max': 0.357}
    extent_0416_dict = {'x_min': -0.044, 'x_max': 0.031,
                        'y_min': -0.022, 'y_max': 0.027,
                        'z_min': 0.416, 'z_max': 0.416}
    extent_dicts = {0.249: extent_0249_dict,
                    0.302: extent_0302_dict,
                    0.357: extent_0357_dict,
                    0.416: extent_0416_dict}
    return extent_dicts


def maximum_all_probe_extent(plane):
    pass
