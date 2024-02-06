
import numpy as np

def debris_calc(impact_toughness, thickness, u):
    '''
    calculates the minimum worst case size of debris that can be withstood by the blade
    (worst case size is minimum as protection will defend against larger debris)
    authored by: Phil Blecher

    Parameters:
    ----------
    impact_toughness : float
        Impact toughness of the material (J/cm^2)
    thickness : float
        Thickness of the blade (mm)
    u : float
        Velocity of the blade (m/s)

    Returns:
    -------
    worst_length : float
        Minimum worst case size of debris that can be withstood by the blade (mm)
    
    
    '''
    max_energy = impact_toughness * thickness/10 # J
    mass = max_energy/0.5*u**2

    volume = mass/1000 # m^3

    worst_length = volume**(1/3) # assuming the debris is a cube as that would be the smallest volume for a given mass

    return worst_length * 1000 # convert to mm

def torque(river_vel, runner_diameter, r_drum, L, RPM):
    '''
    calculates the torque produced on the blade (N*m)
    authored by: Henry Haslam

    Parameters:
    ----------
    river_vel : float
        River velocity (m/s)
    runner_diameter : float
        Diameter of the runner (m)
    r_drum : float
        Radius of the drum (m)
    L : float
        Length of the blade (m)
    RPM : float
        Rotations per minute of the drum

    Returns:
    -------
    T : float
        Torque produced on the blade (N*m)
    
    
    '''
    # Define variables
    H = (runner_diameter/2) - r_drum # Length of blade from drum to tip (m)
    A = L*H*2 # Area of two blades in contact with water (m^2)
    D = r_drum + H/2 # Distance from drum to center of blade (m)

    # Define constants
    rho = 1000 # Density of fresh water (kg/m^3)
    CD = 1.28 # Drag coefficient

    # Calculate relative velocity (at the root of the blade):
    angular_vel = RPM * 2 * np.pi / 60 # Angular velocity (rad/s)
    blade_vel = angular_vel * r_drum # Blade velocity (m/s)
    rel_vel = river_vel - blade_vel # Relative velocity (m/s)

    T = 0.5 * rho * rel_vel**2 * A * CD * D
    return T

def minimum_blade_thickness(sigma_y, river_vel, L, r_drum, rho=1000, RPM=40, CD=1.28):
    '''
    calculates the minimum blade thickness for a given set of parameters
    authored by: Henry Haslam

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