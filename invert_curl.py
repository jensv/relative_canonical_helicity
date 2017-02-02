import numpy as np
import scipy.fftpack as fft


def fourier_inverse_curl(Bx, By, Bz, x, y, z):
    r"""
    Invert curl with pseudo-spectral method described in MacKay 2006.
    """
    shape = Bx.shape

    Bx_padded = np.pad(Bx, pad_width=zip(shape, shape),
                       mode='reflect')
    By_padded = np.pad(By, pad_width=zip(shape, shape),
                       mode='reflect')
    Bz_padded = np.pad(Bz, pad_width=zip(shape, shape),
                       mode='reflect')

    kx1 = np.zeros(Bx_padded[0, :, 0].size)
    ky1 = np.zeros(By_padded[:, 0, 0].size)
    kz1 = np.zeros(Bz_padded[0, 0, :].size)

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

    Bx_k = fft.fftn(Bx_padded)
    By_k = fft.fftn(By_padded)
    Bz_k = fft.fftn(Bz_padded)

    k_squared = kx**2. + ky**2. + kz**2.

    Ax_k = 1j*(ky*Bz_k - kz*By_k)/k_squared
    Ay_k = 1j*(kz*Bx_k - kx*Bz_k)/k_squared
    Az_k = 1j*(kx*By_k - ky*Bx_k)/k_squared

    Ax_k[0, 0, 0] = 0.
    Ay_k[0, 0, 0] = 0.
    Az_k[0, 0, 0] = 0.

    Ax_padded = np.abs(fft.ifftn(Ax_k))
    Ay_padded = np.abs(fft.ifftn(Ay_k))
    Az_padded = np.abs(fft.ifftn(Az_k))

    Ax = Ax_padded[shape[0]:shape[0]*2,
                   shape[1]:shape[1]*2,
                   shape[2]:shape[2]*2]
    Ay = Ay_padded[shape[0]:shape[0]*2,
                   shape[1]:shape[1]*2,
                   shape[2]:shape[2]*2]
    Az = Az_padded[shape[0]:shape[0]*2,
                   shape[1]:shape[1]*2,
                   shape[2]:shape[2]*2]

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
