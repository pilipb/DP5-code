
import numpy as np
import matplotlib.pyplot as plt

def wall_calc(worst_case_v = 3.5, l_wall = 425e-3, h = 165e-3, sigma_fat = 27.25e6, tau_mod = 0.8505e9, rho_mat = 1200, cpkg = 2.865):
    '''
    Calculates the required thickness of the wall, the yield and shear strength of the wall, and the cost of the wall
    authored by: Oliver Legood

    Defaults given are for the material polycarbonate

    Parameters:
    ----------
    worst_case_v : float
        Worst case velocity of the river (m/s)
    l_wall : float
        Length of the wall (m)
    h : float
        Height of the wall at long end (m)
    sigma_fat : float
        Fatigue strength of the material (Pa)
    tau_mod : float
        Shear modulus of the material (Pa)
    rho_mat : float
        Density of the material (kg/m^3)
    cpkg : float
        Cost of the material per kg

    Returns:
    -------
    t_req : float
        Required thickness of the wall (mm)
    yield_strength_mpa : float
        Yield strength of the wall (MPa)
    shear_strength_mpa : float
        Shear strength of the wall (MPa)
    cost_wall : float
        Cost of the wall

    '''

    # Wall Thickness Calcs
    rho = 1000
    CD = 1.28
    w = 0.5 * rho * worst_case_v**2 * CD * l_wall
    M_max = (w * l_wall**2) / 2
    V = w * l_wall

    # Calculate required thickness (in m)
    t_req_fatigue = (12 * M_max) / (sigma_fat * h**2) 
    t_req_shear = (3 * V) / (2 * tau_mod * h) 

    # Round up to nearest mm
    t_req = (max(t_req_fatigue, t_req_shear))
             
    # Calculate mass and cost per wall segment / wall
    A_segment = (425e-3 * 150e-3) + (0.5 * 15e-3 * 425e-3) # tapering section for our final design

    mass_segment = A_segment * t_req * rho_mat
    mass_wall = mass_segment * 6
    cost_wall = mass_wall * cpkg

    I = (t_req * h**3) / 12
    yield_strength = (M_max * h / 2) / I
    yield_strength_mpa = yield_strength / 1e6

    shear_strength = (3 * V) / (2 * h * t_req)
    shear_strength_mpa = shear_strength / 1e6

    return t_req, yield_strength_mpa, shear_strength_mpa, cost_wall





def drum_calc(G=26.5e9, w=3229, l_blade=0.6, l_drum = 1.2, R=200e-3, rho_mat = 2910, cpkg = 3.04):
    '''
    Calculates the shear and yield strength of the drum for a given thickness
    authored by: Oliver Legood

    Defaults given are for the material HDPE

    Parameters:
    ----------
    G : float
        Shear modulus of the material (Pa)
    w : float
        uniform distributed load (N/m)
    l_blade : float
        length of the blade (m)
    l_drum : float
        length of the drum (m)
    R : float
        outer radius of the drum (m)
    rho_mat : float
        density of the material (kg/m^3)
    cpkg : float
        cost of the material per kg

    Returns:
    -------
    thickness : float
        Required thickness of the drum (mm)
    shear_strength_mpa : float
        Shear strength of the drum (MPa)
    yield_strength_mpa : float
        Yield strength of the drum (MPa)
    shear_mod : float
        Shear modulus of the drum (GPa)
    cost : float
        Cost of the drum


    '''
    # Calculate the forces and moments
    F = w * l_blade
    M = F * l_blade / 2
    T = M * l_drum
    thi_max = np.pi / 180
  
    # Finding required thickness in mm
    thickness = (R - (R**4 - ((2 * l_drum * T) / (thi_max * G * np.pi)))**(1 / 4)) * 1e4

    r_in = R - thickness

    # Calculate the shear and yield strength
    J = (np.pi / 2) * (R**4 - r_in**4)
    shear_strength_pa = (T * R) / J
    shear_strength_mpa = shear_strength_pa / 1e6
    shear_mod = ((shear_strength_mpa * l_drum) / (thi_max * R)) / 1e3
    yield_strength_pa = (M * R) / J
    yield_strength_mpa = yield_strength_pa / 1e6

    # cost calcs
    mass = rho_mat * (np.pi * l_drum * (R**2 - (r_in)**2))
    cost = mass * cpkg

    return thickness, shear_strength_mpa, yield_strength_mpa, shear_mod, cost



def froude_number(velocity, depth, g=9.81):
    '''
    calculates the Froude number of the river

    Parameters:
    ----------
    velocity : float
        Velocity of the river (m/s)
    depth : float
        Depth of the river (m)
    g : float
        Acceleration due to gravity (m/s^2)

    Returns:
    -------
    froude : float
        Froude number of the river
    '''
    froude = velocity / np.sqrt(g * depth)
    return froude

def debris_angle(u, l, rho_debris=580):
    '''
    calculates the maximum angle of debris defenders to send debris under the turbine

    Based on:
    0.5 * rho * u^2 * A * CD * cos(theta) > (m * g - rho * g * V ) * sin(theta)
    where:
    rho = density of water (kg/m^3)
    u = velocity of the blade (m/s)
    A = area of the debris in perpendicular to the flow (m^2)
    CD = drag coefficient
    theta = angle of the debris defender (radians)
    m = mass of the debris (kg)
    g = acceleration due to gravity (m/s^2)
    V = volume of the debris (m^3)

    Rearranging:
    tan(theta) < (u^2 * rho * CD) / (2 * g * l * (rho_debris -rho))
    l = V/A = diameter of the prism of debris (e.g. cylinder) (m) (volume/area perpendicular to flow)

    Parameters:
    ----------
    u : float
        Velocity of the river (m/s)
    l : float
        Diameter of the prism of debris (e.g. cylinder) (m)
    rho_debris : float (mean tropical wood density = 580 kg/m^3 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015RG000514 )
        Density of the debris (kg/m^3)

    Returns:
    -------
    theta : float   
        Maximum angle of debris defenders to send debris under the turbine (radians)
    
    '''

    # Define constants
    rho = 1000 # Density of water (kg/m^3)
    c_d = 1 # Drag coefficient of a bluff body
    g = 9.81 # Acceleration due to gravity (m/s^2)

    tan_theta = (u**2 * rho * c_d) / (2 * g * l * (rho_debris - rho))
    theta = np.arctan(tan_theta)
    theta = np.degrees(theta)
    return theta


def pontoon_height_calc(pontoon_area, turbine_mass):
    '''
    calculates the height of the each of two pontoons required to support the turbine
    authored by: Phil Blecher

    Parameters:
    ----------
    pontoon_area : float
        Area of the pontoon (m^2)
    turbine_mass : float
        Mass of the turbine (kg)

    Returns:
    -------
    depth : float
        depth of the pontoon required to support the turbine (m)
    '''
    # Define constants
    rho = 1000 # Density of water (kg/m^3)

    depth = turbine_mass / (pontoon_area * rho)

    return depth/2

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

    rel_vel = calc_rel_vel(river_vel, r_drum, RPM)

    T = 0.5 * rho * rel_vel**2 * A * CD * D
    return T

def power(river_vel, runner_diameter, r_drum, L, RPM):
    '''
    calculates the power produced by the blade (W)

    Parameters:
    ----------
    river_vel : float
        River velocity (m/s)
    runner_diameter : float
        Diameter of the runner (m)
    r_drum : float
        Radius of the drum (m)
    L : float
        width of a single blade (m) (assume half turbine width)
    RPM : float
        Rotations per minute of the drum

    Returns:
    -------
    P : float
        Power produced by the blade (W)
    '''
    # Define variables
    H = (runner_diameter/2) - r_drum # Length of blade from drum to tip (m)
    A = L*H*2 # Area of two blades in contact with water (m^2)
    D = r_drum + H/2 # Distacne from drum to center of blade (m)

    # Define constants
    rho = 1000 # Density of fresh water (kg/m^3)
    CD = 1.28 # Drag coefficient

    # Calculate relative velocity (at the root of the blade):
    angular_vel = RPM * 2 * np.pi / 60 # Angular velocity (rad/s)
    blade_vel = angular_vel * r_drum # Blade velocity (m/s)
    rel_vel = river_vel - blade_vel # Relative velocity (m/s)

    power = 0.5 * rho * rel_vel**2 * A * CD * D * angular_vel

    # Power limiting
    # max_power = 2000 # W
    # try:
    #     power[power > max_power] = max_power
    # except:
    #     if power > max_power:
    #         power = max_power
    
    ############## GET VALUES
    # calculate the loss due to cupping ( assuming that the cupping is a pyramid shape)
    # vol = ((L*runner_diameter/2) * bucket_depth)/3
    # cup_rad = (0.25*runner_diameter/2 + r_drum)
    # cup_loss = vol * cup_rad * rho * angular_vel
    # power -= cup_loss    

    # add the effect of losses due to small holes for the water to escape
    # A = 0.01 # area of the hole per blade(m^2)
    # loss = 0.5 * rho * rel_vel**2 * A * CD * D * angular_vel
    # power -= loss


    return power

def calc_rel_vel(river_vel, r_drum, RPM):
    '''
    calculates the relative velocity of the blade to the water (m/s)

    Parameters:
    ----------
    river_vel : float
        River velocity (m/s)
    r_drum : float
        Radius of the drum (m) or radius of interest
    RPM : float
        Rotations per minute of the drum

    Returns:
    -------
    rel_vel : float
        Relative velocity of the blade to the water (m/s)
    '''
    # Calculate relative velocity (at the root of the blade):
    angular_vel = RPM * 2 * np.pi / 60 # Angular velocity (rad/s)
    blade_vel = angular_vel * r_drum # Blade velocity (m/s)
    rel_vel = river_vel - blade_vel # Relative velocity (m/s)

    return rel_vel

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