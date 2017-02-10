import numpy as np
import scipy.fftpack as fft
import sys
sys.path.append('../laplace_solver/')
import laplace_solver as lsolve
from scipy.integrate import cumtrapz

def fourier_inverse_curl(Bx, By, Bz, x, y, z, method='fourier', pad=True):
    r"""
    Invert curl with pseudo-spectral method described in MacKay 2006.
    """
    shape = Bx.shape
    Bx_copy = np.array(Bx)
    By_copy = np.array(By)
    Bz_copy = np.array(Bz)

    if pad:
        Bx = np.pad(Bx, pad_width=zip(shape, shape),
                    mode='reflect')
        By = np.pad(By, pad_width=zip(shape, shape),
                    mode='reflect')
        Bz = np.pad(Bz, pad_width=zip(shape, shape),
                    mode='reflect')

    kx1 = np.zeros(Bx[0, :, 0].size)
    ky1 = np.zeros(By[:, 0, 0].size)
    kz1 = np.zeros(Bz[0, 0, :].size)

    dx = x[0, 1, 0] - x[0, 0, 0]
    dy = y[1, 0, 0] - y[0, 0, 0]
    dz = z[0, 0, 1] - z[0, 0, 0]

    nx = kx1.size
    ny = ky1.size
    nz = kz1.size

    kx1 = fft.fftfreq(nx, dx)
    ky1 = fft.fftfreq(ny, dy)
    kz1 = fft.fftfreq(nz, dz)

    kx, ky, kz = np.meshgrid(kx1, ky1, kz1)

    if method == 'fourier':
        Bx_k = np.fft.fftn(Bx)
        By_k = np.fft.fftn(By)
        Bz_k = np.fft.fftn(Bz)
    if method == 'cosine':
        Bx_k = lsolve.dct_3d(shape, Bx)
        By_k = lsolve.dct_3d(shape, By)
        Bz_k = lsolve.dct_3d(shape, Bz)

    k_squared = kx**2. + ky**2. + kz**2.

    if method == 'fourier':
        Ax_k = 1j*(ky*Bz_k - kz*By_k)/k_squared
        Ay_k = 1j*(kz*Bx_k - kx*Bz_k)/k_squared
        Az_k = 1j*(kx*By_k - ky*Bx_k)/k_squared
    if method == 'cosine':
        Ax_k = (ky*Bz_k - kz*By_k)/k_squared
        Ay_k = (kz*Bx_k - kx*Bz_k)/k_squared
        Az_k = (kx*By_k - ky*Bx_k)/k_squared

    Ax_k[0, 0, 0] = 0.
    Ay_k[0, 0, 0] = 0.
    Az_k[0, 0, 0] = 0.

    if method == 'fourier':
        Ax = np.real(np.fft.ifftn(Ax_k))
        Ay = np.real(np.fft.ifftn(Ay_k))
        Az = np.real(np.fft.ifftn(Az_k))

    if method == 'cosine':
        Ax = lsolve.idct_3d(shape, Ax_k)
        Ay = lsolve.idct_3d(shape, Ay_k)
        Az = lsolve.idct_3d(shape, Az_k)

    if pad:
        Ax = Ax[shape[0]:shape[0]*2,
                shape[1]:shape[1]*2,
                shape[2]:shape[2]*2]
        Ay = Ay[shape[0]:shape[0]*2,
                shape[1]:shape[1]*2,
                shape[2]:shape[2]*2]
        Az = Az[shape[0]:shape[0]*2,
                shape[1]:shape[1]*2,
                shape[2]:shape[2]*2]

    B0_x = np.mean(Bx_copy)
    B0_y = np.mean(By_copy)
    B0_z = np.mean(Bz_copy)

    A0_x = -(y*B0_z - z*B0_y)/2.
    A0_y = -(z*B0_x - x*B0_z)/2.
    A0_z = -(x*B0_y - y*B0_x)/2.

    Ax = Ax + A0_x
    Ay = Ay + A0_y
    Az = Az + A0_z

    return [Ax, Ay, Az]


def devore_invert_curl(mesh, b_field, with_z=True):
    r"""
    """
    b_field = np.asarray(b_field)
    dz = mesh[2][0, 0, 1] - mesh[2][0, 0, 0]
    z_length = mesh[0].shape[2]
    A_0x, A_0y = devore_A_0(mesh, b_field)
    A_x = np.expand_dims(A_0x, 2)
    A_y = np.expand_dims(A_0y, 2)
    A_x = np.repeat(A_x, z_length, axis=2)
    A_y = np.repeat(A_y, z_length, axis=2)
    A_x += cumtrapz(b_field[1], axis=2, dx=dz, initial=0)
    A_y -= cumtrapz(b_field[0], axis=2, dx=dz, initial=0)
    if with_z:
        A_z = np.zeros(mesh[0].shape)
        return A_x, A_y, A_z
    else:
        return A_x, A_y


def devore_A_0(mesh, b_field):
    r"""
    """
    dx = mesh[0][0, 1, 0] - mesh[0][0, 0, 0]
    dy = mesh[1][1, 0, 0] - mesh[1][0, 0, 0]
    A_0x = -0.5*cumtrapz(b_field[2, :, :, 0], axis=0, dx=dy, initial=0)
    A_0y = 0.5*cumtrapz(b_field[2, :, :, 0], axis=1, dx=dx, initial=0)
    return A_0x, A_0y
