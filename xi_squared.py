import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append('../../vector_calculus')
import vector_calculus as vc


def calc_and_plot_dists(field1, field2,
                        field1_title=None,
                        field2_title=None,
                        units=None):
    r"""
    """
    field1_magnitudes, field1_orientations = sample_vector_field(field1)
    field2_magnitudes, field2_orientations = sample_vector_field(field2)
    field1_histograms = bin_samples(field1_magnitudes, field1_orientations)
    field2_histograms = bin_samples(field2_magnitudes, field2_orientations)
    xi_sq = xi_squared(field1_histograms['counts_2d'], field2_histograms['counts_2d'])
    plot_histograms(field1_magnitudes, field1_orientations, units=units)
    sns.plt.suptitle(field1_title)
    plot_histograms(field2_magnitudes, field2_orientations, units=units)
    sns.plt.suptitle(field2_title)
    plt.show()
    print 'Xi^2 = %f10.2' % xi_sq


def sample_vector_field(vector_field):
    r"""
    """
    vector_field = np.asarray(vector_field)
    orientations = np.ravel(vector_orientation(vector_field))
    magnitudes = np.ravel(vc.magnitude(vector_field))
    return magnitudes, orientations


def bin_samples(magnitudes,
                orientations,
                magnitude_bins=50,
                orientation_bins=50,
                joint_bins=50):
    r"""
    """
    mag_counts, mag_bins = np.histogram(magnitudes,
                                        bins=magnitude_bins)
    o_counts, o_bins = np.histogram(orientations,
                                    bins=orientation_bins)
    (counts_2d, mag_bins_2d,
     o_bins_2d) = np.histogram2d(magnitudes,
                                 orientations,
                                 bins=(magnitude_bins,
                                       orientation_bins))
    histograms = {'mag_counts': mag_counts,
                  'mag_bins': mag_bins,
                  'o_counts': o_counts,
                  'o_bins': o_bins,
                  'counts_2d': counts_2d,
                  'mag_bins_2d': mag_bins_2d,
                  'o_bins_2d': o_bins_2d}
    return histograms


def plot_histograms(magnitudes, orientations, bins=50,
                    color='red', cmap='Reds', units=None):
    r"""
    """
    joint_grid = sns.JointGrid(magnitudes,
                               orientations)
    joint_grid.plot_joint(plt.hist2d, bins=bins,
                          cmap=cmap)
    joint_grid.plot_marginals(sns.distplot,
                              kde=False, bins=bins,
                              color=color)
    xlabel = 'magnitude' + ' [' + units + ']'
    joint_grid.set_axis_labels(xlabel, 'orientation [rad]')
    return joint_grid


def xi_squared(dist1, dist2, dof=2.):
    r"""
    """
    assert dist1.shape == dist2.shape, "Distributions do not have equal dimensions."
    addends = (dist1 -  dist2)**2 / (dist1 + dist2)
    addends[np.isclose(dist1 + dist2, 0)] = 0
    xi_sq =  np.sum(addends)
    if not dof:
        dof = len(dist1.shape)
    xi_sq = xi_sq/dof
    return xi_sq


def vector_orientation(vector_field):
    r"""
    Return vector angle with the z-axis.
    """
    mag = vc.magnitude(vector_field)
    angle = np.arccos(vector_field[2, :, :, :]/mag)
    angle[np.isclose(mag, 0)] = 0
    reflex = np.where(np.logical_and(vector_field[0, :, :, :] < 0,
                                     vector_field[1, :, :, :] >= 0))
    angle[reflex] = 2.*np.pi - angle[reflex]
    reflex = np.where(np.logical_and(vector_field[0, :, :, :] < 0,
                                     vector_field[1, :, :, :] < 0))
    angle[reflex] = 2.*np.pi - angle[reflex]
    return angle
