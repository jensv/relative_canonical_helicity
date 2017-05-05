import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context('poster')

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as axisartist
import matplotlib.patches as patches

from scipy.constants import proton_mass as m_i
from scipy.constants import elementary_charge as q_e
import scipy.fftpack as fft
from scipy.signal import boxcar
from astropy.convolution import convolve


def compare_helicities(helicities,
                       kinetic=True, relative_kinetic=True,
                       cross=True, relative_cross=True,
                       magnetic=True, relative_magnetic=True,
                       S_0_velocity='u_alfven', normalize=None,
                       nondim=False, absolute=False,
                       filter_width_cross=None, filter_width_kinetic=None,
                       add_cross_magnetic=False, add_three=False,
                       ylim=None, labels_case='default', axes=None,
                       scale='linear', patch_ends=True, shift_time=125):

    assert (scale == 'linear' or scale == 'log' or
            scale == 'symlog'), "scale must be one of linear, log or symlog"
    assert (S_0_velocity is None or S_0_velocity == 'u_alfven' or
            S_0_velocity == 'u_sound'), "S_0_velocity must be one of None, u_alfven, u_sound"

    labels = {'kinetic': r'$\mathcal{H}$',
              'relative_kinetic': r'$\mathcal{H}_{rel}$',
              'cross': r'$X$',
              'relative_cross': r'$\mathcal{X}_{rel}$',
              'magnetic': r'$\mathcal{K}$',
              'relative_magnetic': r'$\mathcal{K}_{rel}$'}

    labels_dimless = {'kinetic': r'$\bar{\mathcal{H}}$',
                      'relative_kinetic': r'$\bar{\mathcal{H}}_{rel}$',
                      'cross': r'$\bar{\mathcal{X}}$',
                      'relative_cross': r'$\bar{\mathcal{X}}_{rel}$',
                      'magnetic': r'$\bar{\mathcal{K}}$',
                      'relative_magnetic': r'$\bar{\mathcal{K}}_{rel}$'}


    labels_dict = {'default': labels, 'dimless': labels_dimless}

    labels = labels_dict[labels_case]

    m_0 = 1.67e-27
    q_0 = 1.6e-19
    l_0 = 0.5
    B_0 = 0.02
    u_0_a = 1.4e5
    u_0_cs = 20e3
    n_0 = 1e18


    if S_0_velocity == 'u_alfven':
        S_0 = l_0*q_0*B_0/(m_0*u_0_a)
    elif S_0_velocity == 'u_sound':
        S_0 =  l_0*q_0*B_0/(m_0*u_0_cs)
    else:
        S_0 = 1.


    if nondim:
        kinetic_divider = m_0**2 * u_0_a**2 * l_0**2*n_0**2
        cross_divider = m_0*q_0*u_0_a*B_0*l_0**3*1./S_0*n_0**2
        magnetic_divider = q_0**2*B_0**2*l_0**4*1./S_0**2.*n_0**2
    else:
        kinetic_divider = 1.
        cross_divider = 1.
        magnetic_divider = 1.

    if not axes:
        axes = plt.gca()

    time = helicities['time']

    keys = helicities.keys()
    keys.remove('time')
    helicities_to_plot = {}
    for key in keys:
        if shift_time:
            helicities[key] = np.roll(helicities[key], shift_time)
        if absolute:
            helicities_to_plot[key] = np.abs(helicities[key])
        else:
            helicities_to_plot[key] = helicities[key]



    if normalize:
        norm = np.max(np.abs(helicities_to_plot[normalize]))
        if 'kinetic' in normalize:
            norm /= kinetic_divider
        if 'cross' in normalize:
            norm /= cross_divider
        if 'magnetic' in normalize:
            norm /= magnetic_divider
    else:
        norm = 1


    if kinetic:
        axes.plot(time, helicities_to_plot['kinetic']/(norm*kinetic_divider),
                c='blue', ls='--', label=labels['kinetic'])
    if relative_kinetic:
        if filter_width_kinetic:
            box = boxcar(filter_width_kinetic)/filter_width_kinetic
            signal = np.asarray(helicities_to_plot['relative_kinetic']/(norm*kinetic_divider))
            signal_filt = convolve(signal, box, boundary='wrap')
            kinetic_final = np.asarray(signal_filt)
            axes.plot(time, signal_filt, c='blue', ls='-', label=labels['relative_kinetic'])
            axes.plot(time, signal,
                      c='blue', alpha=0.2)
        else:
            kinetic_final = helicities_to_plot['relative_kinetic']/(norm*kinetic_divider)
            axes.plot(time, helicities_to_plot['relative_kinetic']/(norm*kinetic_divider),
                      c='blue', ls='-', label=labels['relative_kinetic'])
    if cross:
        axes.plot(time, helicities_to_plot['cross']/(norm*cross_divider),
                 c='green', ls='--', label=labels['cross'])
    if relative_cross:
        if filter_width_cross:
            box = boxcar(filter_width_cross)/filter_width_cross
            signal = np.asarray(helicities_to_plot['relative_cross']/(norm*cross_divider))
            signal_filt = convolve(signal, box, boundary='wrap')
            cross_final = np.asarray(signal_filt)
            axes.plot(time, signal_filt, c='green', ls='-', label=labels['relative_cross'])
            axes.plot(time, signal,
                      c='green', alpha=0.4)
        else:
            cross_final = helicities_to_plot['relative_cross']/(norm*cross_divider)
            axes.plot(time, helicities_to_plot['relative_cross']/(norm*cross_divider),
                      c='green', ls='-', label=labels['relative_cross'])
    if magnetic:
        axes.plot(time, helicities_to_plot['magnetic']/(norm*magnetic_divider),
                 c='red', ls='--', label=labels['magnetic'])
    if relative_magnetic:
        magnetic_final = helicities_to_plot['relative_magnetic']/(norm*magnetic_divider)
        axes.plot(time, helicities_to_plot['relative_magnetic']/(norm*magnetic_divider),
                 c='red', ls='-', label=labels['relative_magnetic'])
    axes.set_xlabel(r'$t$ [$\mu s$]')
    axes.set_yscale(scale)
    axes.set_ylabel(r'$K$ [$J$ $kg$ $m^2$]')
    if normalize:
        axes.set_ylabel(r'$K$ [-]')

    if ylim:
        axes.set_ylim(ylim)
    if add_cross_magnetic:
        axes.plot(time, cross_final + magnetic_final,
                  c='yellow', ls='-', label=labels['relative_magnetic'] + " $+$ " + labels['relative_cross'])
    if add_three:
        axes.plot(time, kinetic_final + cross_final + magnetic_final,
                  c='black', ls='-', label=labels['relative_magnetic'] + " $+$ " + labels['relative_cross'] + " $+$ " + labels['relative_kinetic'])
    axes.legend(loc='best', fancybox=True, frameon=True, framealpha=0.9)
    if patch_ends:
        in_dark_box_1 = patches.Rectangle((5.644, -1000),
                                          11.9-5.644, 2000., alpha=0.4, color='grey')
        in_light_box_1 = patches.Rectangle((0.748, -1000),
                                           2.108-0.748, 2000., alpha=0.1, color='grey')
        in_light_box_2 = patches.Rectangle((2.584, -1000),
                                           12.104-2.584, 2000, alpha=0.1, color='grey')
        axes.add_patch(in_dark_box_1)
        axes.add_patch(in_light_box_1)
        axes.add_patch(in_light_box_2)
    return axes


def compare_helicities_mean_std(helicities_mean, helicities_std,
                                kinetic=True, relative_kinetic=True,
                                cross=True, relative_cross=True,
                                magnetic=True, relative_magnetic=True,
                                S_0_velocity='u_alfven', normalize=None,
                                nondim=False, absolute=False,
                                filter_width_cross=None, filter_width_kinetic=None,
                                add_cross_magnetic=False, add_three=False,
                                ylim=None, labels_case='default', axes=None,
                                scale='linear', patch_ends=True, shift_time=125):
    r"""
    """
    assert (scale == 'linear' or scale == 'log' or
            scale == 'symlog'), "scale must be one of linear, log or symlog"
    assert (S_0_velocity is None or S_0_velocity == 'u_alfven' or
            S_0_velocity == 'u_sound'), "S_0_velocity must be one of None, u_alfven, u_sound"

    labels = {'kinetic': r'$\mathcal{H}$',
              'relative_kinetic': r'$\mathcal{H}_{rel}$',
              'cross': r'$X$',
              'relative_cross': r'$\mathcal{X}_{rel}$',
              'magnetic': r'$\mathcal{K}$',
              'relative_magnetic': r'$\mathcal{K}_{rel}$'}

    labels_dimless = {'kinetic': r'$\bar{\mathcal{H}}$',
                      'relative_kinetic': r'$\bar{\mathcal{H}}_{rel}$',
                      'cross': r'$\bar{\mathcal{X}}$',
                      'relative_cross': r'$\bar{\mathcal{X}}_{rel}$',
                      'magnetic': r'$\bar{\mathcal{K}}$',
                      'relative_magnetic': r'$\bar{\mathcal{K}}_{rel}$'}


    labels_dict = {'default': labels, 'dimless': labels_dimless}

    labels = labels_dict[labels_case]

    m_0 = 1.67e-27
    q_0 = 1.6e-19
    l_0 = 0.5
    B_0 = 0.02
    u_0_a = 1.4e5
    u_0_cs = 20e3
    n_0 = 1e18


    if S_0_velocity == 'u_alfven':
        S_0 = l_0*q_0*B_0/(m_0*u_0_a)
    elif S_0_velocity == 'u_sound':
        S_0 =  l_0*q_0*B_0/(m_0*u_0_cs)
    else:
        S_0 = 1.


    if nondim:
        kinetic_divider = m_0**2 * u_0_a**2 * l_0**2*n_0**2
        cross_divider = m_0*q_0*u_0_a*B_0*l_0**3*1./S_0*n_0**2
        magnetic_divider = q_0**2*B_0**2*l_0**4*1./S_0**2.*n_0**2
    else:
        kinetic_divider = 1.
        cross_divider = 1.
        magnetic_divider = 1.

    if not axes:
        axes = plt.gca()

    time = helicities_mean['time']

    keys = helicities_mean.keys()
    keys.remove('time')
    helicities_to_plot = {}
    for key in keys:
        if shift_time:
            helicities_mean[key] = np.roll(helicities_mean[key], shift_time)
            helicities_std[key] = np.roll(helicities_std[key], shift_time)
        if absolute:
            helicities_to_plot[key] = np.abs(helicities_mean[key])
        else:
            helicities_to_plot[key] = helicities_mean[key]



    if normalize:
        norm = np.max(np.abs(helicities_to_plot[normalize]))
        if 'kinetic' in normalize:
            norm /= kinetic_divider
        if 'cross' in normalize:
            norm /= cross_divider
        if 'magnetic' in normalize:
            norm /= magnetic_divider
    else:
        norm = 1


    if kinetic:
        axes.plot(time, helicities_to_plot['kinetic']/(norm*kinetic_divider),
                  c='blue', ls='--', label=labels['kinetic'])
        axes.fill_between(time,
                          (helicities_to_plot['kinetic'] - helicities_std['kinetic'])/(norm*kinetic_divider),
                          (helicities_to_plot['kinetic'] + helicities_std['kinetic'])/(norm*kinetic_divider),
                          facecolor='#95d0fc')
    if relative_kinetic:
        if filter_width_kinetic:
            box = boxcar(filter_width_kinetic)/filter_width_kinetic
            signal = np.asarray(helicities_to_plot['relative_kinetic']/(norm*kinetic_divider))
            signal_filt = convolve(signal, box, boundary='wrap')
            kinetic_final = np.asarray(signal_filt)
            axes.plot(time, signal_filt, c='blue', ls='-', label=labels['relative_kinetic'])
            axes.plot(time, signal,
                      c='blue', alpha=0.2)
        else:
            kinetic_final = helicities_to_plot['relative_kinetic']/(norm*kinetic_divider)
            kinetic_std_final = helicities_std['relative_kinetic']/(norm*kinetic_divider)
            axes.plot(time, helicities_to_plot['relative_kinetic']/(norm*kinetic_divider),
                      c='blue', ls='-', label=labels['relative_kinetic'])
            axes.fill_between(time,
                              (helicities_to_plot['relative_kinetic'] - helicities_std['relative_kinetic'])/(norm*kinetic_divider),
                              (helicities_to_plot['relative_kinetic'] + helicities_std['relative_kinetic'])/(norm*kinetic_divider),
                              facecolor='#95d0fc')

    if cross:
        axes.plot(time, helicities_to_plot['cross']/(norm*cross_divider),
                 c='green', ls='--', label=labels['cross'])
        axes.fill_between(time,
                          (helicities_to_plot['cross'] - helicities_std['cross'])/(norm*cross_divider),
                          (helicities_to_plot['cross'] + helicities_std['cross'])/(norm*cross_divider),
                          facecolor='#2dfe54')
    if relative_cross:
        if filter_width_cross:
            box = boxcar(filter_width_cross)/filter_width_cross
            signal = np.asarray(helicities_to_plot['relative_cross']/(norm*cross_divider))
            signal_filt = convolve(signal, box, boundary='wrap')
            cross_final = np.asarray(signal_filt)
            axes.plot(time, signal_filt, c='green', ls='-', label=labels['relative_cross'])
            axes.plot(time, signal,
                      c='green', alpha=0.4)
        else:
            cross_final = helicities_to_plot['relative_cross']/(norm*cross_divider)
            cross_std_final = helicities_std['relative_cross']/(norm*cross_divider)
            axes.plot(time, helicities_to_plot['relative_cross']/(norm*cross_divider),
                      c='green', ls='-', label=labels['relative_cross'])
            axes.fill_between(time,
                              (helicities_to_plot['relative_cross'] - helicities_std['relative_cross'])/(norm*cross_divider),
                              (helicities_to_plot['relative_cross'] + helicities_std['relative_cross'])/(norm*cross_divider),
                              facecolor='#2dfe54')
    if magnetic:
        axes.plot(time, helicities_to_plot['magnetic']/(norm*magnetic_divider),
                 c='red', ls='--', label=labels['magnetic'])
        axes.fill_between(time,
                          (helicities_to_plot['magnetic'] - helicities_std['magnetic'])/(norm*magnetic_divider),
                          (helicities_to_plot['magnetic'] + helicities_std['magnetic'])/(norm*magnetic_divider),
                          facecolor='#ff474c')
    if relative_magnetic:
        magnetic_final = helicities_to_plot['relative_magnetic']/(norm*magnetic_divider)
        magnetic_std_final = helicities_std['relative_magnetic']/(norm*magnetic_divider)
        axes.plot(time, helicities_to_plot['relative_magnetic']/(norm*magnetic_divider),
                 c='red', ls='-', label=labels['relative_magnetic'])
        axes.fill_between(time,
                          (helicities_to_plot['relative_magnetic'] - helicities_std['relative_magnetic'])/(norm*magnetic_divider),
                          (helicities_to_plot['relative_magnetic'] + helicities_std['relative_magnetic'])/(norm*magnetic_divider),
                          facecolor='#ff474c')
    axes.set_xlabel(r'$t$ [$\mu s$]')
    axes.set_yscale(scale)
    axes.set_ylabel(r'$K$ [$J$ $kg$ $m^2$]')
    if normalize:
        axes.set_ylabel(r'$K$ [-]')

    if ylim:
        axes.set_ylim(ylim)
    if add_cross_magnetic:
        axes.plot(time, cross_final + magnetic_final,
                  c='yellow', ls='-', label=labels['relative_magnetic'] + " $+$ " + labels['relative_cross'])
    if add_three:
        std_propagated = np.sqrt(kinetic_final**2*kinetic_std_final**2 +
                                 cross_final**2*cross_std_final**2 +
                                 magnetic_final**2*magnetic_std_final**2)
        total_helicity = kinetic_final + cross_final + magnetic_final
        axes.plot(time, kinetic_final + cross_final + magnetic_final,
                  c='black', ls='-', label=labels['relative_magnetic'] + " $+$ " + labels['relative_cross'] + " $+$ " + labels['relative_kinetic'])
        axes.fill_between(time,
                          (total_helicity - std_propagated),
                          (total_helicity + std_propagated),
                          facecolor='black', alpha=0.5)


    axes.legend(loc='best', fancybox=True, frameon=True, framealpha=0.9)
    if patch_ends:
        in_dark_box_1 = patches.Rectangle((5.644, -1000), 11.9-5.644, 2000.,
                                          alpha=0.4, color='grey')
        in_light_box_1 = patches.Rectangle((0.748, -1000), 2.108-0.748, 2000.,
                                           alpha=0.1, color='grey')
        in_light_box_2 = patches.Rectangle((2.584, -1000), 12.104-2.584, 2000,
                                           alpha=0.1, color='grey')
        axes.add_patch(in_dark_box_1)
        axes.add_patch(in_light_box_1)
        axes.add_patch(in_light_box_2)
    return axes
