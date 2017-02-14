import numpy as np
from invert_curl.invert_curl import devore_invert_curl
from laplace_solver.laplace_solver import laplace_3d_dct_fd
import vector_calculus.vector_calculus as vc


def determine_reference_fields(mesh, circulation):
    r"""
    Return reference fields used for relative helicity.
    Reference fields consist of a circulation and the general momentum vector
    of which the curl gives the circulation.

    Parameters
    ----------
    mesh: list of ndarray
        3D mesh

    circulation: list of ndarray
        3D vector which is the curl quantity.
        e.g. magnetic field B or flow vorticity omega.

    Returns
    -------
    momentum_ref: list of ndarray
        reference general momentum field e.g.
        reference magnetic vector potential or
        reference flow
    circulation_ref: list of ndarray
        curl of reference field e.g.
        reference magnetic field, reference flow vorticity.
    Notes
    -----
    Circulation reference dotted with surface normal should be
    the negative of the real circulation dotted with the surface
    normal.

    .. math::
    $-\vec{Circ}_{ref} \cdot \hat{n}= \vec{Circ} \cdot \hat{n}$
    """
    boundary = make_boundary(circulation)
    scalar_potential_ref = laplace_3d_dct_fd(mesh, boundary)
    circulation_ref = vc.gradient(scalar_potential_ref, mesh=mesh)
    momentum_ref = devore_invert_curl(mesh,
                                      circulation_ref)
    return momentum_ref, circulation_ref


def make_boundary(field):
    r"""
    Return boundary conditions for circulation reference.

   .. math::
    $-\vec{Circ}_{ref} \cdot \hat{n}= \vec{Circ} \cdot \hat{n}$
    """
    boundary = np.zeros(field[0].shape)
    boundary[:, 0, :] = -field[0][:, 0, :]
    boundary[:, -1, :] = -field[0][:, -1, :]
    boundary[0, :, :] = -field[1][0, :, :]
    boundary[-1, :, :] = -field[1][-1, :, :]
    boundary[:, :, 0] = -field[2][:, :, 0]
    boundary[:, :, -1] = -field[2][:, :, -1]
    return boundary
