
import numpy as np

def minimum_blade_thickness(sigma_y, river_vel, L, r_drum, rho=1000, RPM=40, CD=1.28):
    '''
    calculates the minimum blade thickness for a given set of parameters

    Parameters:
    ----------
    sigma_y : float
        Yield strength of the material (Pa)
    river_vel : float
        River velocity (m/s)
    L : float
        Length of the blade (m)
    r_drum : float
        Radius of the drum (m)
    rho : float
        Density of water (kg/m^3)
    RPM : float
        Rotations per minute of the drum
    CD : float
        Drag coefficient of the blade

    Returns:
    -------
    h : float
        Minimum blade thickness (mm)
    '''
    # Calculate relative velocity (at the root of the blade):
    angular_vel = RPM * 2 * np.pi / 60 # Angular velocity (rad/s)
    blade_vel = angular_vel * r_drum # Blade velocity (m/s)
    rel_vel = river_vel - blade_vel # Relative velocity (m/s)

    # Calculate minimum blade thickness (mm):
    h = ((rho * (rel_vel**2)*CD*L**2)/(4 * sigma_y))**0.5*1000

    return h