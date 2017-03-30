import numpy as np
import reference_field as ref
from scipy.integrate import trapz
from scipy.constants import proton_mass as m_i
from scipy.constants import elementary_charge as q_e


factor = {'kinetic': m_i**2,
          'cross': q_e*m_i,
          'magnetic': q_e**2}

def kinetic_helicity(velocity,
                     vorticity,
                     dx, dy, dz,
                     density=None):
    r"""
    Returns kinetic helicity.
    """
    return general_helicity(velocity,
                            vorticity,
                            'kinetic',
                            dx, dy, dz,
                            density=density)


def cross_helicity(velocity,
                   magnetic_field,
                   dx, dy, dz, density=None):
    r"""
    Returns cross helicity.
    """
    return 2.*general_helicity(velocity,
                               magnetic_field,
                               'cross',
                               dx, dy, dz,
                               density=density)


def magnetic_helicity(vector_potential,
                      magnetic_field,
                      dx, dy, dz, density=None):
    r"""
    Returns magnetic helicity.
    """
    return general_helicity(vector_potential,
                            magnetic_field,
                            'magnetic',
                            dx, dy, dz,
                            density=density)


def general_helicity(momentum, circulation, kind,
                     dx, dy, dz, density=None):
    r"""
    Returns a general helicity for a given momentum and circulation.
    """
    momentum = np.asarray(momentum)
    ciruclation = np.asarray(circulation)
    helicity_density = np.einsum('ijkl, ijkl -> jkl',
                                 momentum, circulation)
    helicity_density = helicity_density*factor[kind]
    helicity = volume_integral_trapz(helicity_density,
                                     dx, dy, dz,
                                     density=density)
    return helicity


def rel_kin_helicity(velocity, velocity_ref,
                         vorticity, vorticity_ref,
                         dx, dy, dz, density=None):
    r"""
    Returns relative kinetic helicity.
    """
    return general_relative_helicity(velocity, velocity_ref,
                                     vorticity, vorticity_ref,
                                     'kinetic', dx, dy, dz,
                                     density=density)


def rel_cross_helicity(velocity, velocity_ref,
                       magnetic_field, magnetic_field_ref,
                       dx, dy, dz, density=None):
    r"""
    Returns relative cross helicity.

    Notes
    -----

    Dot product differs from other relative helicities.
    u_- \cdot B_+ + u_+ \cdot B_-
    """
    term1 = general_relative_helicity(velocity, velocity_ref,
                                      magnetic_field, magnetic_field_ref,
                                      'cross', dx, dy, dz,
                                      density=density)
    inverted_velocity_ref = [-velocity_ref[0], -velocity_ref[1],
                             -velocity_ref[2]]
    inverted_magnetic_field_ref = [-magnetic_field_ref[0], -magnetic_field_ref[1],
                                   -magnetic_field_ref[1]]
    term2 = general_relative_helicity(velocity, inverted_velocity_ref,
                                      magnetic_field, inverted_magnetic_field_ref,
                                      'cross', dx, dy, dz,
                                      density=density)
    return term1 + term2


def rel_mag_helicity(vector_potential, vector_potential_ref,
                     magnetic_field, magnetic_field_ref,
                     dx, dy, dz, density=None):
    r"""
    Returns relative magnetic helicity.
    """
    return general_relative_helicity(vector_potential, vector_potential_ref,
                                     magnetic_field, magnetic_field_ref,
                                     'magnetic', dx, dy, dz,
                                     density=density)



def general_relative_helicity(momentum, momentum_ref,
                              circulation, circulation_ref,
                              kind,
                              dx, dy, dz, density=None):
    r"""
    Returns general relative helicity.

    Notes
    -----
    Dot product is
    (momentum - momentum_ref) \cdot (circulation + circulation_ref)
    """
    momentum = np.asarray(momentum)
    momentum_ref = np.asarray(momentum_ref)
    circulation = np.asarray(circulation)
    ciruclation_ref = np.asarray(circulation_ref)
    relative_momentum = momentum - momentum_ref
    relative_circulation = circulation + circulation_ref
    return general_helicity(relative_momentum, relative_circulation,
                            kind, dx, dy, dz,
                            density=density)


def volume_integral_trapz(quantity, dx, dy, dz, density=None):
    r"""

    Notes
    -----
    Each trapezoidal integration reduces the number of dimensions by 1.
    """
    if density is None:
        density = 1.
    integral = trapz(trapz(trapz(quantity*density**2,
                                 dx=dy, axis=0),
                           dx=dx, axis=0),
                     dx=dz, axis=0)
    return integral
