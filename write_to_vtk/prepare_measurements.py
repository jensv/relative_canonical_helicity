import numpy as np
import scipy.io.idl as idl
from mach_probe_analysis import ion_current_to_mach_number as ic_to_mach
from read_from_sql import read_from_sql

def read_idl(quantity, data_path='../output/comprehensive_3d_plot/2016-08-12/'):
    r"""
    Read idl files for all planes.
    """
    idl_ending = '.sav'
    z249_file = data_path + quantity + '_z249' + idl_ending
    z302_file = data_path + quantity + '_z302' + idl_ending
    z357_file = data_path + quantity + '_z357' + idl_ending
    z416_file = data_path + quantity + '_z416' + idl_ending

    z249_measurements = idl.readsav(z249_file)
    z249_measurements['a_out'] = z249_measurements['a_out'].astype('float64')
    z302_measurements = idl.readsav(z302_file)
    z302_measurements['a_out'] = z302_measurements['a_out'].astype('float64')
    z357_measurements = idl.readsav(z357_file)
    z357_measurements['a_out'] = z357_measurements['a_out'].astype('float64')
    z416_measurements = idl.readsav(z416_file)
    z416_measurements['a_out'] = z416_measurements['a_out'].astype('float64')

    measurements = {0.249: z249_measurements,
                    0.302: z302_measurements,
                    0.357: z357_measurements,
                    0.416: z416_measurements}

    return measurements


def remove_points_outside_box(measurements,
                              x_min= -0.025, x_max=0.025,
                              y_min=-0.016, y_max=0.017,
                              z_min=0.249, z_max=0.416):
    r"""
    Return measurements dict with points outside cube of given bounds removed.
    """
    points = np.dstack((measurements['x_out'],
                        measurements['y_out'],
                        measurements['z_out']))[0]
    outside_box = np.where(np.logical_or.reduce((points[:, 0] < x_min,
                                                    points[:, 0] > x_max,
                                                    points[:, 1] < y_min,
                                                    points[:, 1] > y_max,
                                                    points[:, 2] < z_min,
                                                    points[:, 2] > z_max)))[0]
    for key in ['x_out', 'y_out', 'z_out']:
        measurements[key] = np.delete(measurements[key],
                                      outside_box)
    for key in ['std', 'a_out']:
        for time_point in xrange(measurements['delays'].size):
            measurements[key][time_point] = np.delete(measurements[key][time_point],
                                                      outside_box)
    return measurements


def average_duplicate_points(data_dict):
    r"""
    Find duplicate points, average them and record standard deviation.
    """
    data_dict['x_out'] = data_dict['x_out'].astype('float64')
    data_dict['y_out'] = data_dict['y_out'].astype('float64')
    data_dict['a_out'] = data_dict['a_out'].astype('float64')
    time_points = data_dict['a_out'].shape[0]
    data = {}
    for idx in xrange(data_dict['x_out'].size):
        location = (data_dict['x_out'][idx], data_dict['y_out'][idx])
        if location in data.keys():
            data[location] = np.column_stack((data[location], data_dict['a_out'][:, idx]))
        else:
            data[location] = data_dict['a_out'][:, idx]

    unique_data_dict = {'x_out': [],
                        'y_out': [],
                        'a_out': [],
                        'std': []}
    for location in data.keys():
        if data[location][0].size > 1:
            unique_data_dict['std'].append(data[location].std(axis=1, ddof=1))
            unique_data_dict['a_out'].append(data[location].mean(axis=1))
        else:
            unique_data_dict['std'].append(np.zeros(time_points))
            unique_data_dict['a_out'].append(data[location])
        unique_data_dict['x_out'].append(location[0])
        unique_data_dict['y_out'].append(location[1])

    unique_data_dict['x_out'] = np.asarray(unique_data_dict['x_out'])
    unique_data_dict['y_out'] = np.asarray(unique_data_dict['y_out'])
    unique_data_dict['a_out'] = np.hsplit(np.asarray(unique_data_dict['a_out']), time_points)
    unique_data_dict['std'] = np.hsplit(np.asarray(unique_data_dict['std']), time_points)
    unique_data_dict['delays'] = data_dict['delays']
    return unique_data_dict


def read_points_from_measurement_dict(measurement_dict, time_point, z_planes):
    r"""
    Collectand return measurement points and values from dictionary.
    """
    x_points = np.empty((0))
    y_points = np.empty((0))
    z_points = np.empty((0))
    values = np.empty((0))
    for z_plane in z_planes:
        plane_measurements = measurement_dict[z_plane]
        x_points = np.append(x_points, plane_measurements['x_out'])
        y_points = np.append(y_points, plane_measurements['y_out'])
        z_points = np.append(z_points, np.ones(plane_measurements['x_out'].size)*z_plane)
        values = np.append(values, plane_measurements['a_out'][time_point])

    points = [x_points, y_points, z_points]
    points = np.asarray(points)
    points = np.swapaxes(points, 0, 1)
    return points, values


def read_mach_probe_data(args):
    r"""
    """
    timesteps = args.mach_time_steps
    database = args.shot_database
    table = args.table_name
    min_spectral_density = args.min_spectral

    z_direction_1, z_direction_2 = 0, 180
    y_direction_1, y_direction_2 = 90, 270
    angle_signs = {0: 1,
                   180: -1,
                   90: -1,
                   0: 1}

    condition_z_0416 = ("campaigns = 'mach_probe_plane_campaign_1'"
                        " AND fiducial_pre_crowbar_gyration_spectral_density > "
                        + str(min_spectral_density) +
                        " AND mach_signals_exist = 1"
                        " AND (mach_orientation = " + str(z_direction_1) +
                        " OR mach_orientation = " + str(z_direction_2) + ")")

    condition_y_0416 = ("campaigns = 'mach_probe_plane_campaign_1'"
                        " AND fiducial_pre_crowbar_gyration_spectral_density > "
                        + str(min_spectral_density) +
                        " AND mach_signals_exist = 1"
                        " AND (mach_orientation = " + str(y_direction_1) +
                        " OR mach_orientation = " + str(y_direction_2) + ")")

    cursor, connection = read_from_sql.cursor_with_rows(condition_z_0416,
                                                        database,
                                                        table)
    z_0416_shots = cursor.fetchall()
    cursor.close()
    connection.close()

    cursor, connection = read_from_sql.cursor_with_rows(condition_y_0416,
                                                        database,
                                                        table)
    y_0416_shots = cursor.fetchall()
    cursor.close()
    connection.close()

    condition_z_302 = ("campaigns = 'mach_probe_plane_campaign_2'"
                       " AND fiducial_pre_crowbar_gyration_spectral_density > "
                       + str(min_spectral_density) +
                       " AND mach_signals_exist = 1"
                       " AND (mach_orientation = " + str(z_direction_1) +
                       " OR mach_orientation = " + str(z_direction_2) + ")")

    cursor, connection = read_from_sql.cursor_with_rows(condition_z_302,
                                                        database,
                                                        table)
    z_0302_shots = cursor.fetchall()
    cursor.close()
    connection.close()

    mach_z_0416_measurements = ic_to_mach.run_mach_analysis(z_0416_shots,
                                                            timesteps,
                                                            angle_signs)
    mach_y_0416_measurements = ic_to_mach.run_mach_analysis(y_0416_shots,
                                                            timesteps,
                                                            angle_signs)
    mach_z_0302_measurements = ic_to_mach.run_mach_analysis(z_0302_shots,
                                                            timesteps,
                                                            angle_signs)

    mach_z_0416_measurements['delays'] = np.arange(timesteps)
    mach_y_0416_measurements['delays'] = np.arange(timesteps)
    mach_z_0302_measurements['delays'] = np.arange(timesteps)

    mach_y_measurements = {0.416: mach_y_0416_measurements}
    mach_z_measurements = {0.416: mach_z_0416_measurements}

    return mach_y_measurements, mach_z_measurements


def cut_and_average_quantity(measurements, box_extent, planes, bounds=None):
    for plane in planes:
        measurements[plane] = average_duplicate_points(measurements[plane])
    if bounds:
        measurements = remove_points_out_of_bounds(measurements,
                                                    bounds[0],
                                                    bounds[1],
                                                    planes)
    all_planes = combine_all_planes(measurements, planes)
    if box_extent:
       all_planes = remove_points_outside_box(all_planes, box_extent[0],
                                              box_extent[1], box_extent[2],
                                              box_extent[3], box_extent[4],
                                              box_extent[5])
    all_planes = remove_nan_points(all_planes)
    return all_planes


def remove_points_out_of_bounds(data_dict, lower, upper, planes):
    r"""
    Find points out of bounds and remove them.
    """
    new_data_dict = {}
    for plane in planes:
        new_data_dict[plane] = {'x_out': [],
                                'y_out': [],
                                'a_out': [],
                                'std': []}
        to_remove = []
        for time in xrange(len(data_dict[plane]['a_out'])):
            indexes = (np.where(np.logical_or(data_dict[plane]['a_out'][time] < lower,
                                              data_dict[plane]['a_out'][time] > upper))[0])
            if indexes.size > 0:
                to_remove.append(indexes)
        if len(to_remove) > 0:
            to_remove = np.concatenate(to_remove)
        to_remove = np.unique(to_remove)
        new_data_dict[plane]['x_out'] = np.delete(data_dict[plane]['x_out'], to_remove)
        new_data_dict[plane]['y_out'] = np.delete(data_dict[plane]['y_out'], to_remove)
        for time in xrange(len(data_dict[plane]['a_out'])):
            new_data_dict[plane]['a_out'].append(np.delete(data_dict[plane]['a_out'][time],
                                                 to_remove))
            new_data_dict[plane]['std'].append(np.delete(data_dict[plane]['std'][time],
                                                         to_remove))
        new_data_dict[plane]['a_out'] = np.asarray(new_data_dict[plane]['a_out'])
        new_data_dict[plane]['std'] = np.asarray(new_data_dict[plane]['std'])
        new_data_dict[plane]['delays'] = data_dict[plane]['delays']
    return new_data_dict


def combine_all_planes(measurements, planes):
    r"""
    """
    all_planes = dict(measurements[planes[0]])
    all_planes['z_out'] = planes[0]*np.ones(all_planes['x_out'].size)
    for plane in planes[1:]:
        for key in ['x_out', 'y_out']:
            all_planes[key] = np.concatenate((all_planes[key],
                                              measurements[plane][key]))
        all_planes['z_out'] = np.concatenate((all_planes['z_out'],
                                              plane*np.ones(measurements[plane]['x_out'].size)))
        for key in ['std', 'a_out']:
            for time_point in xrange(measurements[plane]['delays'].size):
                all_planes[key][time_point] = np.concatenate((all_planes[key][time_point],
                                                              measurements[plane][key][time_point]))
    return all_planes


def remove_nan_points(measurements):
    r"""
    """
    nan_positions = np.where(np.isnan(measurements['a_out']))[1]
    for key in ['x_out', 'y_out', 'z_out']:
        measurements[key] = np.delete(measurements[key],
                                      nan_positions)
    for key in ['std', 'a_out']:
        measurements[key] = np.delete(measurements[key],
                                      nan_positions, axis=1)
    return measurements


def remove_positions(measurements, positions_to_remove):
    r"""
    """
    for key in ['x_out', 'y_out', 'z_out']:
        measurements[key] = np.delete(measurements[key],
                                      positions_to_remove)
    for key in ['std', 'a_out']:
        measurements[key] = np.delete(measurements[key],
                                      positions_to_remove, axis=1)
    return measurements


def remove_plane(plane, measurements):
    plane_positions = np.where(measurements['z_out'] == plane)[0]
    for key in ['x_out', 'y_out', 'z_out']:
        measurements[key] = np.delete(measurements[key],
                                      plane_positions)
    for key in ['std', 'a_out']:
        measurements[key] = np.delete(measurements[key],
                                      plane_positions, axis=1)
    return measurements
