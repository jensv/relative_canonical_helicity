import numpy as np
import scipy.sparse as sparse
import scipy.linalg as linalg
import scipy.fftpack as fft

def laplace_3d_dct_fd(mesh, boundary_values):
    r"""
    """
    d_x = mesh[0][0, 1, 0] - mesh[0][0, 0, 0]
    d_y = mesh[1][1, 0, 0] - mesh[1][0, 0, 0]
    d_z = mesh[2][0, 0, 1] - mesh[2][0, 0, 0]
    deltas = [d_x, d_y, d_z]
    shape = mesh[0].shape
    right_side = construct_right_side(shape, deltas, boundary_values)
    transform_right_side = dct_3d(shape, right_side)
    eigenvalues = construct_eigenvalues(shape)
    transform_solution = transform_right_side/eigenvalues
    solution = idct_3d(shape, transform_solution)
    return solution


def construct_right_side(shape, deltas, boundary_values):
    r"""
    """
    right_side = np.zeros(shape)
    right_side[:, 0, :] = 2.*deltas[0]*boundary_values[:, 0, :]
    right_side[0, :, :] = 2.*deltas[1]*boundary_values[0, :, :]
    right_side[:, :, 0] = 2.*deltas[2]*boundary_values[:, :, 0]
    right_side[:, -1, :] = -2.*deltas[0]*boundary_values[:, -1, :]
    right_side[-1, :, :] = -2.*deltas[1]*boundary_values[-1, :, :]
    right_side[:, :, -1] = -2.*deltas[2]*boundary_values[:, :, -1]
    return right_side


def dct_3d(shape, scalar):
    r"""
    """
    transform = fft.dct(fft.dct(fft.dct(scalar, axis=1, type=1),
                                axis=0, type=1),
                        axis=2, type=1)
    return transform


def idct_3d(shape, transform):
    r"""
    """
    ifft_factor_x = 0.5/(shape[1] - 1.)
    ifft_factor_y = 0.5/(shape[0] - 1.)
    ifft_factor_z = 0.5/(shape[2] - 1.)
    scalar = ifft_factor_x*fft.dct(transform, axis=1, type=1)
    scalar = ifft_factor_y*fft.dct(scalar, axis=0, type=1)
    scalar = ifft_factor_z*fft.dct(scalar, axis=2, type=1)
    return scalar


def construct_eigenvalues(shape):
    r"""
    """
    x_index = np.arange(0, shape[1])
    y_index = np.arange(0, shape[0])
    z_index = np.arange(0, shape[2])
    index_mesh = np.meshgrid(x_index,
                             y_index,
                             z_index)
    eigenvalues = 2.*(np.cos(np.pi*index_mesh[0] /
                             (shape[1] - 1.)) +
                      np.cos(np.pi*index_mesh[1] /
                             (shape[0] - 1.)) +
                      np.cos(np.pi*index_mesh[2] /
                             (shape[2] - 1.)) - 3.)
    eigenvalues[np.isclose(eigenvalues, 0, atol=1e-12)] = 1.
    return eigenvalues
