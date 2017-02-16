# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 13:48:25 2015

@author: Jens von der Linden
"""

import numpy as np
from os.path import exists
import scipy.io.idl as idl
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
import visit_writer

def save_to_unstructured_grid(measurements,
                              visit_quantity_name,
                              out_dir,
                              prefix='_unstructured_grid_'):
    r"""
    """
    file_prefix = out_dir + '/' + visit_quantity_name + prefix
    assert (len(measurements['x_out']) ==
            len(measurements['y_out']) ==
            len(measurements['z_out']))
    (points,
     connectivity,
     variables_all_time) = prepare_for_unstructured_vtk(measurements,
                                                        visit_quantity_name)
    #assert len(points) == len(measurements['x_out'])*3
    #assert len(variables_all_time[0][0][3]) == len(measurements['x_out'])
    write_all_time_unstructured(file_prefix, points,
                                connectivity, variables_all_time)

def prepare_for_unstructured_vtk(measurements, quantity_name):
    r"""
    """
    if np.unique(np.asarray(measurements['z_out'])).size > 1:
        points = np.dstack((measurements['x_out'],
                            measurements['y_out'],
                            measurements['z_out']))[0]
        dimensions = 3
    else:
        points = np.dstack((measurements['x_out'],
                            measurements['y_out']))[0]
        dimensions = 2
    triangulation = Delaunay(points)
    if dimensions == 2:
        points = np.dstack((measurements['x_out'],
                            measurements['y_out'],
                            0.419*np.ones(len(measurements['x_out']))))[0]
    points = tuple(points.ravel())
    if dimensions == 3:
        connectivity = tuple([(visit_writer.tetrahedron, int(simplex[0]),
                               int(simplex[1]), int(simplex[2]), int(simplex[3]))
                               for simplex in triangulation.simplices])
    else:
        connectivity = tuple([(visit_writer.triangle, int(simplex[0]),
                               int(simplex[1]), int(simplex[2]))
                               for simplex in triangulation.simplices])
    variables_all_time = []
    for time_point in xrange(measurements['delays'].size):
        variables = ((quantity_name, 1, 1, tuple(measurements['a_out'][time_point])),)
        variables_all_time.append(variables)
    return points, connectivity, variables_all_time

def  write_all_time_unstructured(file_prefix, points, connectivity, variables_all_time):
     r"""
     """
     for time_point in xrange(len(variables_all_time)):
         path = file_prefix + str(time_point).zfill(4)
         visit_writer.WriteUnstructuredMesh(path, 1, points,
                                            connectivity,
                                            variables_all_time[time_point])
