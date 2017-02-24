import numpy as np
from scipy.special import ellipe, ellipk

def curl(quantity, spacing=(1, 1, 1), mesh=None,
         vector_grad=None):
    r"""
    Return 3D curl.
    """
    if not vector_grad:
        dx, dy, dz = spacing
        if mesh:
            dx = mesh[0][0, 1, 0] - mesh[0][0, 0, 0]
            dy = mesh[1][1, 0, 0] - mesh[1][0, 0, 0]
            dz = mesh[2][0, 0, 1] - mesh[2][0, 0, 0]
        dx_dy = np.gradient(quantity[0], axis=0)/dy
        dx_dz = np.gradient(quantity[0], axis=2)/dz
        dy_dx = np.gradient(quantity[1], axis=1)/dx
        dy_dz = np.gradient(quantity[1], axis=2)/dz
        dz_dx = np.gradient(quantity[2], axis=1)/dx
        dz_dy = np.gradient(quantity[2], axis=0)/dy
    else:
        dx_dy = vector_grad[0][1]
        dx_dz = vector_grad[0][2]
        dy_dx = vector_grad[1][0]
        dy_dz = vector_grad[1][2]
        dz_dx = vector_grad[2][0]
        dz_dy = vector_grad[2][1]
    curl_x = dz_dy - dy_dz
    curl_y = dx_dz - dz_dx
    curl_z = dy_dx - dx_dy
    return curl_x, curl_y, curl_z


def dot_product(vector1, vector2):
    r"""
    Return dot product.
    """
    assert vector1[0].shape == vector2[0].shape, 'Vector fields do not have the same dimensions.'
    assert 4 > len(vector1) > 1, 'Vectors should have at least 2 no more then 3 components.'
    if len(vector1) == 3:
        return vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2]
    else:
        return vector1[0]*vector2[0] + vector1[1]*vector2[1]


def gradient(scalar, dx=1, dy=1, dz=1,
             mesh=None):
    if mesh:
        dx = mesh[0][0, 1, 0] - mesh[0][0, 0, 0]
        dy = mesh[1][1, 0, 0] - mesh[1][0, 0, 0]
        dz = mesh[2][0, 0, 1] - mesh[2][0, 0, 0]
    grad = np.gradient(scalar, dy, dx, dz)
    return grad[1], grad[0], grad[2]


def magnitude(vector):
    r"""
    Return magnitude of vector.
    """
    assert 4 > len(vector) > 1, 'Vectors should have at least 2 no more then 3 components.'
    if len(vector) == 3:
        return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    else:
        return np.sqrt(vector[0]**2 + vector[1]**2)


def field_from_wire(limits=(-5, 5, -5, 5, -5, 5),
                    points=(11, 11, 11),
                    center=(0, 0),
                    mu_0I=1, a=1, wire_along_axis='y'):
    r"""
    """
    x, y, z = np.meshgrid(np.linspace(limits[0], limits[1], points[0]),
                          np.linspace(limits[2], limits[3], points[1]),
                          np.linspace(limits[4], limits[5], points[2]))
    mesh = [x, y, z]
    B_x = np.zeros(x.shape)
    B_y = np.zeros(x.shape)
    B_z = np.zeros(x.shape)
    A_x = np.zeros(x.shape)
    A_y = np.zeros(x.shape)
    A_z = np.zeros(x.shape)

    if wire_along_axis == 'y':
        theta = np.arctan2(z - center[1], x - center[0])
        r = np.sqrt((x - center[0])**2. + (z - center[1])**2.)
    elif wire_along_axis == 'z':
        theta = np.arctan2(y - center[1], x - center[0])
        r = np.sqrt((x - center[0])**2. + (y - center[1])**2.)

    inside = np.where(r <= a)
    outside = np.where(r > a)
    inside = (inside[0], inside[1], inside[2])
    outside = (outside[0], outside[1], outside[2])

    if wire_along_axis == 'y':
        B_x[inside] = -mu_0I*r[inside]/(2.*np.pi*a**2)*np.sin(theta[inside])
        B_z[inside] = mu_0I*r[inside]/(2.*np.pi*a**2)*np.cos(theta[inside])
        B_x[outside] = -mu_0I/(2.*np.pi*r[outside])*np.sin(theta[outside])
        B_z[outside] = mu_0I/(2.*np.pi*r[outside])*np.cos(theta[outside])

        A_y[inside] = -mu_0I/(4.*np.pi*a**2)*(r[inside]**2 - a**2)
        A_y[outside] = -mu_0I/(2.*np.pi)*np.log(r[outside]/a)
    elif wire_along_axis == 'z':
        B_x[inside] = -mu_0I*r[inside]/(2.*np.pi*a**2)*np.sin(theta[inside])
        B_y[inside] = mu_0I*r[inside]/(2.*np.pi*a**2)*np.cos(theta[inside])
        B_x[outside] = -mu_0I/(2.*np.pi*r[outside])*np.sin(theta[outside])
        B_y[outside] = mu_0I/(2.*np.pi*r[outside])*np.cos(theta[outside])

        A_z[inside] = -mu_0I/(4.*np.pi*a**2)*(r[inside]**2 - a**2)
        A_z[outside] = -mu_0I/(2.*np.pi)*np.log(r[outside]/a)

    return mesh, A_x, A_y, A_z, B_x, B_y, B_z


def field_from_loop(limits=(-5, 5, -5, 5, -5, 5),
                    points=(11, 11, 11),
                    mu_0I=1., a=1.):
    r"""
    """
    x, y, z = np.meshgrid(np.linspace(limits[0], limits[1], points[0]),
                          np.linspace(limits[2], limits[3], points[1]),
                          np.linspace(limits[4], limits[5], points[2]))
    mesh = [x, y, z]
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)

    B_x = np.zeros(x.shape)
    B_y = np.zeros(x.shape)
    B_z = np.zeros(x.shape)
    A_x = np.zeros(x.shape)
    A_y = np.zeros(x.shape)
    A_z = np.zeros(x.shape)

    k_sq = (4.*a*r)/((a+r)**2+z**2)
    E = ellipe(k_sq)
    K = ellipk(k_sq)
    A_phi = mu_0I/(np.pi*np.sqrt(k_sq))*np.sqrt(a/r)*((1 - k_sq/2.)*K - E)
    B_z = (mu_0I/(2*np.pi)*1./np.sqrt((a + r)**2 + z**2)*
           (K + (a**2 - r**2 - z**2)/((a - r)**2 + z**2)*E))
    B_r = (mu_0I/(2*np.pi)*z/(r*np.sqrt((a + r)**2 + z**2))*
           (-K + (a**2 + r**2 + z**2)/((a - r)**2 + z**2)*E))



    B_z =  (1./np.sqrt((a + r)**2 + z**2)*
            (K + E*(a**2 - r**2 - z**2)/((a - r)**2 + z**2)))

    B_r = (z/(r*np.sqrt((a + r)**2 + z**2))*
           (-K + E*(a**2 + r**2 + z**2)/((a - r)**2 + z**2)))

    A_x = -A_phi*np.sin(phi)
    A_y = A_phi*np.cos(phi)
    B_x = B_r*np.cos(phi)
    B_y = B_r*np.sin(phi)

    return mesh, A_x, A_y, A_z, B_x, B_y, B_z


def constant_B_field(limits=(-5, 5, -5, 5, -5, 5),
                     points=(11, 11, 11),
                     direction='y'):
    r"""
    """
    x, y, z = np.meshgrid(np.linspace(limits[0], limits[1], points[0]),
                          np.linspace(limits[2], limits[3], points[1]),
                          np.linspace(limits[4], limits[5], points[2]))
    mesh = [x, y, z]

    B_x = np.zeros(x.shape)
    B_y = np.zeros(x.shape)
    B_z = np.zeros(x.shape)

    if direction == 'x':
        B_x = np.ones(x.shape)
    elif direction == 'y':
        B_y = np.ones(x.shape)
    elif direction == 'z':
        B_z = np.ones(x.shape)

    return mesh, B_x, B_y, B_z

