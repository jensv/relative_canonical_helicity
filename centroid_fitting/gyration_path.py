r""""
Created on Apr 17 2017

@author: Jens von der Linden

Plot field line null gyration of RSX magnetic flux tube.
"""

import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_style('white')
sns.set_context('poster')

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as axisartist
import matplotlib.patches as patches

def gyration_path(axes=None, circles=None, step=25,
                  start=0, error_test=False,
                  points=None, errors=None,
                  field_null_path="/home/jensv/rsx/jens_analysis/output/field_nulls/",
                  field_null_file="2017-04-13-23-41/field_nulls.txt",
                  measurement_limits=(-0.022, 0.024, -0.017, 0.018),
                  bxby_limits = (-0.024, 0.025, -0.073, 0.038),
                  errorevery=10, params_guess=[0, 0, 0.01],
                  circle_fit=False, xlim=(-0.03, 0.05),
                  xticks=None):

    centroid_file = field_null_path + field_null_file
    if not points is None:
        field_nulls = points
    else:
        field_nulls = np.loadtxt(centroid_file)

    x_min, x_max = measurement_limits[0], measurement_limits[1]
    y_min, y_max = measurement_limits[2], measurement_limits[3]

    if not axes:
        fig, axes = plt.subplots(1, 1)

    measurement_box = patches.Rectangle((x_min, y_min), x_max-x_min, y_max-y_min,
                                    color='grey', alpha=0.4)
    bx_by_x_min, bx_by_x_max = bxby_limits[0:2]
    bx_by_y_min, bx_by_y_max = bxby_limits[2:4]
    bx_by_measurement_box = patches.Rectangle((bx_by_x_min, bx_by_y_min),
                                              bx_by_x_max - bx_by_x_min,
                                              bx_by_y_max - bx_by_y_min,
                                              color='grey', alpha=0.1)
    axes.add_patch(measurement_box)
    axes.add_patch(bx_by_measurement_box)

    colormap = np.linspace(0, 1, 250-start)
    axes.scatter(field_nulls[start:, 0], field_nulls[start:, 1], c=colormap, zorder=100)
    if not errors is None:
        axes.errorbar(field_nulls[start:, 0], field_nulls[start:, 1],
                      xerr=errors[start:, 0], yerr=errors[start:, 1],
                      ecolor='grey',
                      fmt='none', zorder=0, errorevery=errorevery, alpha=0.5)


    #axes.text(-0.008, -0.015, r'$0 \mu s$')
    #axes.text(0.03, -0.003, r'$%2.1f \mu s$' % (0.068*56))
    #axes.text(-0.03, 0.017, r'$%2.1f \mu s$' % (0.068*208))

    if circles:
        for i, field_null in enumerate(field_nulls[::step]):
            colormap = np.linspace(1, 0, np.round(250./step))
            circle = patches.Circle(field_null, radius=0.02, facecolor='none',
                                    edgecolor=str(colormap[i]), alpha=0.5)
            axes.scatter(field_null[0], field_null[1], c='red')
            axes.add_patch(circle)

    if circle_fit:
        circle_params, success = leastsq(to_min, params_guess,
                                         args=np.asarray([field_nulls[:, 0],
                                                          field_nulls[:, 1]]))
        circle = patches.Circle((circle_params[0],  circle_params[1]),
                                radius=circle_params[2], facecolor='none',
                                edgecolor='black', lw=5)
        axes.add_patch(circle)


    axes.set_xlabel('x [m]')
    axes.set_ylabel('y [m]')
    axes.set_xlim(xlim)
    if not xticks is None:
        axes.xaxis.set_ticks(xticks)
    axes.set_aspect('equal')
    axes.invert_xaxis()
    return axes

def to_min(params, points):
    r"""
    Returns circle expression to minimize with least squares.
    """
    a = 2.*params[0]
    b = 2.*params[1]
    c = params[2]**2 - params[1]**2 - params[0]**2
    return a*points[0] + b*points[1] + c - points[0]**2 - points[1]**2
