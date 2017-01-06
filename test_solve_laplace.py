import numpy as np
from laplace_solver import laplace_3d_dct_fd

def test_laplace_3d_dct_fd_solver():
    r"""
    """
    mesh, A_field, B_field = field_from_wire()
    d_x = mesh[0][0, 1, 0] - mesh[0][0, 0, 0]
    d_y = mesh[1][1, 0, 0] - mesh[1][0, 0, 0]
    d_z = mesh[2][0, 0, 1] - mesh[2][0, 0, 0]
    boundary_values = np.zeros(mesh[0].shape)
    boundary_values[:, 0, :] = B_field[0][:, 0, :]
    boundary_values[0, :, :] = B_field[1][0, :, :]
    boundary_values[:, :, 0] = B_field[2][:, :, 0]
    boundary_values[:, -1, :] = B_field[0][:, -1, :]
    boundary_values[-1, :, :] = B_field[1][-1, :, :]
    boundary_values[:, :, -1] = B_field[2][:, :, -1]
    solution = laplace_3d_dct_fd(mesh, boundary_values)
    B_x_from_solution = np.gradient(solution, axis=1)/d_x
    B_y_from_solution = np.gradient(solution, axis=0)/d_y
    B_z_from_solution = np.gradient(solution, axis=2)/d_z
    size = mesh[0].size
    error = np.sum(np.sqrt((B_field[0] - B_x_from_solution)**2.))/size
    error += np.sum(np.sqrt((B_field[1] - B_y_from_solution)**2.))/size
    error += np.sum(np.sqrt((B_field[2] - B_z_from_solution)**2.))/size
    error /= 3.
    print error
    assert error < 1e-3, ("Error in solving known case of "
                          "current through wire larger than expected.")

def field_from_wire():
    r"""
    """
    x = np.linspace(-5, 5, 30)
    y = np.linspace(-5, 5, 20)
    z = np.linspace(-5, 5, 10)

    mesh = np.meshgrid(x, y, z)
    delta_x = x[1] - x[0]
    delta_y = y[1] - y[0]
    delta_z = z[1] - z[0]
    rho_x = 1./delta_x
    rho_y = 1./delta_y
    rho_z = 1./delta_z

    shape = mesh[0].shape

    r = np.sqrt((mesh[0]-15)**2. + mesh[1]**2.)
    current = 1.
    wire_radius = 1.
    mu_0 = 1.
    A_x = np.zeros(mesh[0].shape)
    A_y = np.zeros(mesh[0].shape)
    A_z = -current*mu_0/(2.*np.pi)*np.log(r/wire_radius)
    A_field = [A_x, A_y, A_z]

    theta = np.arctan2(mesh[1], mesh[0] - 15.)
    B_x = mu_0*current/(2.*np.pi*r)*np.sin(theta)*-1
    B_y = mu_0*current/(2.*np.pi*r)*np.cos(theta)
    B_z = np.zeros(mesh[0].shape)
    B_field = [B_x, B_y, B_z]

    return mesh, A_field, B_field
