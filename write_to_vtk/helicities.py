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
                     dx, dy, dz):
    r"""
    """
    return general_helicity(velocity,
                            vorticity,
                            'kinetic',
                            dx, dy, dz)


def cross_helicity(velocity,
                   magnetic_field,
                   dx, dy, dz):
    r"""
    """
    return general_helicity(velocity,
                            magnetic_field,
                            'cross',
                            dx, dy, dz)


def magnetic_helicity(vector_potential,
                      magnetic_field,
                      dx, dy, dz):
    r"""
    """
    return general_helicity(vector_potential,
                            magnetic_field,
                            'magnetic',
                            dx, dy, dz)


def general_helicity(momentum, circulation, kind,
                     dx, dy, dz):
    r"""
    """
    momentum = np.asarray(momentum)
    ciruclation = np.asarray(circulation)
    helicity_density = np.einsum('ijkl, ijkl -> jkl',
                                 momentum, circulation)
    helicity_density = helicity_density*factor[kind]
    helicity = volume_integral_trapz(helicity_density,
                                     dx, dy, dz)
    return helicity


def rel_kin_helicity(velocity, velocity_ref,
                         vorticity, vorticity_ref,
                         dx, dy, dz):
    r"""
    """ 
    return general_relative_helicity(velocity, velocity_ref,
                                     vorticity, vorticity_ref,
                                     'kinetic', dx, dy, dz)


def rel_cross_helicity(velocity, velocity_ref,
                            magnetic_field, magnetic_field_ref,
                            dx, dy, dz):
    r"""
    """
    return general_relative_helicity(velocity, velocity_ref,
                                     magnetic_field, magnetic_field_ref,
                                     'cross', dx, dy, dz)


def rel_mag_helicity(vector_potential, vector_potential_ref,
                               magnetic_field, magnetic_field_ref,
                               dx, dy, dz):
    r"""
    """
    return general_relative_helicity(vector_potential, vector_potential_ref,
                                     magnetic_field, magnetic_field_ref,
                                     'magnetic', dx, dy, dz)



def general_relative_helicity(momentum, momentum_ref,
                              circulation, circulation_ref,
                              kind,
                              dx, dy, dz):
    r"""
    """
    momentum = np.asarray(momentum)
    momentum_ref = np.asarray(momentum_ref)
    circulation = np.asarray(circulation)
    ciruclation_ref = np.asarray(circulation_ref)
    relative_momentum = momentum - momentum_ref
    relative_circulation = circulation + circulation_ref
    return general_helicity(relative_momentum, relative_circulation,
                            kind, dx, dy, dz)


def volume_integral_trapz(quantity, dx, dy, dz):
    r"""

    Notes
    -----
    Each trapezoidal integration reduces the number of dimensions by 1.
    """
    integral = trapz(trapz(trapz(quantity,
                                 dx=dy, axis=0),
                           dx=dx, axis=0),
                     dx=dz, axis=0)
    return integral