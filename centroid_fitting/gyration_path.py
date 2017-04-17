r""""
Created on Apr 17 2017

@author: Jens von der Linden

Plot field line null gyration of RSX magnetic flux tube.
"""

import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_style('white')
sns.set_context('poster')

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as axisartist
import matplotlib.patches as patches

def gyration_path(axes=None, circles=None, step=25,
                  start=0, error_test=False,
                  field_null_path="/home/jensv/rsx/jens_analysis/output/field_nulls/",
                  field_null_file="2017-04-13-23-41/field_nulls.txt",
                  measurement_limits=(-0.022, 0.024, -0.017, 0.018),
                  bxby_limits = (-0.024, 0.025, -0.073, 0.038)):

    centroid_file = field_null_path + field_null_file
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
    if error_test:
        axes.errorbar(field_nulls[start:, 0], field_nulls[start:, 1],
                      xerr=0.01, yerr=0.01, ecolor='grey',
                      fmt='none', zorder=0, errorevery=10, alpha=0.5)


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
    axes.set_xlabel('x [m]')
    axes.set_ylabel('y [m]')
    axes.set_xlim(-0.03, 0.04)
    axes.set_aspect('equal')
    axes.invert_xaxis()
    return axes
