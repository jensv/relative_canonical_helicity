# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 22:30:08 2015

@author: jensv
"""

import numpy as np
import scipy.fftpack as fftpack
from scipy.optimize import curve_fit

from datetime import datetime
import os

import matplotlib.pylab as plt
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context('poster')


def fit_signal(time, signal, pad=0):
    r"""
    Fit signal to a cosine function using initial guesses from fourier analysis
    and min/max finding.

    Parameters
    ----------
    time : ndarray of floats
        sample times
    signal : ndarray of floats
        samples
    pad : int, optional
        number of samples to pad to interpolate between fourier bins
 
    Returns
    -------
    fit_parameters : dict
        dictionary of cosine fit parameters: amplitude, phase_shift, offset, frequency 
    fit_stds : dict
        dictionary of uncertainty in fit parameters.
    covariance : ndarray of float (MxM)
        covariance matrix of fit

    Example
    -------
    >>> time = np.linspace(0., 1./(50e3)*10., 2000) 
    >>> signal = 1.*np.cos(50e3*time*2.*np.pi + 0.5) + 1.
    >>> (fit_parameters, fit_stds, 
         covariance) = fit_signal(time, signal)
    >>> print fit_parameters
        {}
    """
    amplitude_guess = estimate_amplitude(signal)
    offset_guess = estimate_offset(signal, amplitude_guess)
    frequency_guess = estimate_frequency(signal, time, pad=pad)
    phase_guess = estimate_phase(signal, time, frequency_guess)
    (fit_parameters, fit_stds,
     covariance) = cosine_fit(time, signal,
                              amplitude_guess=amplitude_guess,
                              frequency_guess=frequency_guess,
                              offset_guess=offset_guess,
                              phase_guess=phase_guess)
    return fit_parameters, fit_stds, covariance


def fourier_transform(signal, time, pad=0):
    r"""
    Return shifted fourier tansform and frequencies.
    Optionally pad the input signal with zeros.
    """
    padded_signal = np.pad(signal, pad, mode='reflect')
    fourier = fftpack.fft(padded_signal)
    freqs = fftpack.fftfreq(padded_signal.size, d=time[1] - time[0])
    fourier_shifted = fftpack.fftshift(fourier)
    freqs_shifted = fftpack.fftshift(freqs)
    return fourier_shifted, freqs_shifted


def determine_transform_params(fourier, freqs):
    r"""
    Return parameters that can be determined from fourier transform.
    power, amplitude, phase and DC offset.
    """
    n_over_two = freqs.size / 2
    max_index = np.abs(fourier[n_over_two+1:]).argmax()
    frequency = freqs[n_over_two+1:][max_index]
    power = np.abs(fourier[n_over_two+1:][max_index])
    amplitude = power / n_over_two
    phase = np.angle(fourier[n_over_two+1:][max_index])
    DC_offset = np.abs(fourier[n_over_two])/freqs.size
    return frequency, power, amplitude, phase, DC_offset


def estimate_frequency(signal, time, pad=0):
    r"""
    Estimate frequency from maximum of Fourier transform.
    """
    transform, freqs = fourier_transform(signal, time, pad=pad)
    frequency_estimate = determine_transform_params(transform, freqs)[0]
    return frequency_estimate


def estimate_amplitude(signal):
    r"""
    Estiamte amplitude by taking the difference of the min and max of the
    signal.
    """
    maximum = signal.max()
    minimum = signal.min()
    amplitude = (maximum - minimum)/2.
    return amplitude


def estimate_offset(signal, amplitude):
    r"""
    Estiamte offset by taking the difference of the amplitude, max and min.
    """
    offset_from_max = signal.max() - amplitude
    offset_from_min = signal.min() + amplitude
    return np.mean([offset_from_max, offset_from_min])


def estimate_phase(signal, time, frequency):
    r"""output_path
    Estiamte phase of the signal by looking for a max (phase=0) and determing
    how much of a period back the signal start is.
    """
    zero_phase_index = signal.argmax()
    delta_time = time[zero_phase_index] - time[0]
    period = 1./frequency
    period_fraction = 1. - (delta_time % period / period)
    return period_fraction*2.*np.pi


def cosine_to_fit(t, amplitude, frequency, offset, phase):
    r"""
    Cosine function with fit parameters.
    """
    return amplitude*np.cos(t*frequency*2.*np.pi + phase) + offset


def cosine_fit(time, signal, amplitude_guess=1., frequency_guess=1.,
               offset_guess=1., phase_guess=1.):
    r"""
    Least-squares fit cosine function to the signal.
    """
    guess = [amplitude_guess, frequency_guess, offset_guess, phase_guess]
    optimum, covariance = curve_fit(cosine_to_fit, time, signal, guess)
    fit_params = {'amplitude': optimum[0], 'frequency': optimum[1],
                  'offset': optimum[2], 'phase': optimum[3]}
    stds = np.sqrt(np.diag(covariance))
    fit_std = {'amplitude': stds[0], 'frequency': stds[1],
               'offset': stds[2], 'phase': stds[3]}
    return fit_params, fit_std, covariance


def plot_signal_vs_fit(time, signal, fit_params, number):
    r"""
    Plot fit vs signal for visual verification.
    """
    date = datetime.today().strftime('%Y-%m-%d')
    output_path = '../output/' + date
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    fig = plt.figure()
    axes = fig.gca()
    axes.plot(time, signal, label='signal')
    fit = (fit_params['amplitude'] * np.cos(fit_params['frequency']*2.*np.pi*
                                            time + fit_params['phase']) +
           fit_params['offset'])
    axes.plot(time, fit, label='fit')
    axes.legend(loc='best')
    axes.set_xlabel('time [s]')
    axes.set_ylabel('Amplitude')
    fig.savefig(output_path + '/' + number + '.png')
