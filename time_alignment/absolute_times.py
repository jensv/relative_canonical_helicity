# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 10:05:08 2016

@author: Jens von der Linden
"""

import numpy as np


def absolute_times(shot, shot_settings, delays, reference='phase',
                   number_of_delays=None):
    r"""
    Return absolute times for a shot, look up reference in sql database.

    Parameters
    ----------
    shot : int
        shot number
    database_path : string
        path to sql database
    delays : interable
        relative times if phase in radian, if ramp elaped time since start of  ramp
    reference : string
        the reference time to use either phase, or ramp
    number_of_delays : int
        optional - if phase is used instead of supplying specific delays
        this option will return an number of equally spaced times across one
        period.
    Returns
    --------
    times : ndarray of float
        absolute times which correpond to relative times requested.
    """
    msg = "reference has to be 'phase' or 'ramp'"
    assert reference == 'phase' or reference == 'ramp', msg

    if reference == 'phase':
        if number_of_delays:
            delays = np.arange(number_of_delays)*2*np.pi/number_of_delays
        zero_phase_time, period = (shot_settings['zero_phase_time'],
                                   shot_settings['period'])
        times = times_from_phase(zero_phase_time, period, delays)
    else:
        ramp_time = shot_settings['ramp_time']
        times = times_from_current_ramp(ramp_time, delays)

    return times


def times_from_phase(reference_time, period, delays):
    r"""
    """
    delays = np.asarray(delays)
    time_radian = period / (2. * np.pi)
    times = delays*time_radian + reference_time
    return times


def times_from_current_ramp(reference_time, delays):
    r"""
    """
    delays = np.asarray(delays)
    times = reference_time + delays
    return times
