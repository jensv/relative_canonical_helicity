import numpy as np
import scipy.fftpack as fft
import sys
sys.path.append('../laplace_solver/')
import laplace_solver as lsolve

def fourier_inverse_curl(Bx, By, Bz, x, y, z, method='fourier'):
    r"""
    Invert curl with pseudo-spectral method described in MacKay 2006.
    """
    shape = Bx.shape

    if method == 'fourier':
        Bx_padded = np.pad(Bx, pad_width=zip(shape, shape),
                           mode='reflect')
        By_padded = np.pad(By, pad_width=zip(shape, shape),
                           mode='reflect')
        Bz_padded = np.pad(Bz, pad_width=zip(shape, shape),
                           mode='reflect')

        kx1 = np.zeros(Bx_padded[0, :, 0].size)
        ky1 = np.zeros(By_padded[:, 0, 0].size)
        kz1 = np.zeros(Bz_padded[0, 0, :].size)

    if method == 'cosine':
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
        Bx_k = fft.fftn(Bx_padded)
        By_k = fft.fftn(By_padded)
        Bz_k = fft.fftn(Bz_padded)
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
        Ax_padded = np.real(fft.ifftn(Ax_k))
        Ay_padded = np.real(fft.ifftn(Ay_k))
        Az_padded = np.real(fft.ifftn(Az_k))

        Ax = Ax_padded[shape[0]:shape[0]*2,
                       shape[1]:shape[1]*2,
                       shape[2]:shape[2]*2]
        Ay = Ay_padded[shape[0]:shape[0]*2,
                       shape[1]:shape[1]*2,
                       shape[2]:shape[2]*2]
        Az = Az_padded[shape[0]:shape[0]*2,
                       shape[1]:shape[1]*2,
                       shape[2]:shape[2]*2]
    if method == 'cosine':
        Ax = lsolve.idct_3d(shape, Ax_k)
        Ay = lsolve.idct_3d(shape, Ay_k)
        Az = lsolve.idct_3d(shape, Az_k)

    B0_x = np.mean(Bx)
    B0_y = np.mean(By)
    B0_z = np.mean(Bz)

    A0_x = -(y*B0_z - z*B0_y)/2.
    A0_y = -(z*B0_x - x*B0_z)/2.
    A0_z = -(x*B0_y - y*B0_x)/2.

    Ax = Ax + A0_x
    Ay = Ay + A0_y
    Az = Az + A0_z

    return [Ax, Ay, Az]
