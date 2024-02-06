# import libraries and modules needed
import os
import numpy 
import numpy as np
import math
from scipy import integrate, linalg
from matplotlib import pyplot as plt

from aerofoil import naca_foil

class Panel:
    """
    Contains information related to a panel.
    """
    def __init__(self, xa, ya, xb, yb):
        """
        Initializes the panel.
        
        Sets the end-points and calculates the center-point, length,
        and angle (with the x-axis) of the panel.
        Defines if the panel is located on the upper or lower surface of the geometry.
        Initializes the source-strength, tangential velocity, and pressure coefficient
        of the panel to zero.
        
        Parameters
        ---------_
        xa: float
            x-coordinate of the first end-point.
        ya: float
            y-coordinate of the first end-point.
        xb: float
            x-coordinate of the second end-point.
        yb: float
            y-coordinate of the second end-point.
        """
        self.xa, self.ya = xa, ya  # panel starting-point
        self.xb, self.yb = xb, yb  # panel ending-point
        
        self.xc, self.yc = (xa + xb) / 2, (ya + yb) / 2  # panel center
        self.length = numpy.sqrt((xb - xa)**2 + (yb - ya)**2)  # panel length
        
        # orientation of panel (angle between x-axis and panel's normal)
        if xb - xa <= 0.0:
            self.beta = numpy.arccos((yb - ya) / self.length)
        elif xb - xa > 0.0:
            self.beta = numpy.pi + numpy.arccos(-(yb - ya) / self.length)
        
        # panel location
        if self.beta <= numpy.pi:
            self.loc = 'upper'  # upper surface
        else:
            self.loc = 'lower'  # lower surface
        
        self.sigma = 0.0  # source strength
        self.vt = 0.0  # tangential velocity
        self.cp = 0.0  # pressure coefficient

class Freestream:
    """
    Freestream conditions.
    """
    def __init__(self, u_inf=1.0, alpha=0.0):
        """
        Sets the freestream speed and angle (in degrees).
        
        Parameters
        ----------
        u_inf: float, optional
            Freestream speed;
            default: 1.0.
        alpha: float, optional
            Angle of attack in degrees;
            default 0.0.
        """
        self.u_inf = u_inf
        self.alpha = numpy.radians(alpha)  # degrees to radians

def compute_tangential_velocity(panels, freestream, gamma, A_source, B_vortex):
    """
    Computes the tangential surface velocity.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
    freestream: Freestream object
        Freestream conditions.
    gamma: float
        Circulation density.
    A_source: 2D Numpy array of floats
        Source contribution matrix for the normal velocity.
    B_vortex: 2D Numpy array of floats
        Vortex contribution matrix for the normal velocity.
    """
    A = numpy.empty((panels.size, panels.size + 1), dtype=float)
    # matrix of source contribution on tangential velocity
    # is the same than
    # matrix of vortex contribution on normal velocity
    A[:, :-1] = B_vortex
    # matrix of vortex contribution on tangential velocity
    # is the opposite of
    # matrix of source contribution on normal velocity
    A[:, -1] = -numpy.sum(A_source, axis=1)
    # freestream contribution
    b = freestream.u_inf * numpy.sin([freestream.alpha - panel.beta 
                                      for panel in panels])
    
    strengths = numpy.append([panel.sigma for panel in panels], gamma)
    
    tangential_velocities = numpy.dot(A, strengths) + b
    
    for i, panel in enumerate(panels):
        panel.vt = tangential_velocities[i]

def compute_pressure_coefficient(panels, freestream):
    """
    Computes the surface pressure coefficients.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
    freestream: Freestream object
        Freestream conditions.
    """
    for panel in panels:
        panel.cp = 1.0 - (panel.vt / freestream.u_inf)**2

def build_freestream_rhs(panels, freestream):
    """
    Builds the right-hand side of the system 
    arising from the freestream contribution.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
    freestream: Freestream object
        Freestream conditions.
    
    Returns
    -------
    b: 1D Numpy array of floats
        Freestream contribution on each panel and on the Kutta condition.
    """
    b = numpy.empty(panels.size + 1, dtype=float)
    # freestream contribution on each panel
    for i, panel in enumerate(panels):
        b[i] = -freestream.u_inf * numpy.cos(freestream.alpha - panel.beta)
    # freestream contribution on the Kutta condition
    b[-1] = -freestream.u_inf * (numpy.sin(freestream.alpha - panels[0].beta) +
                                 numpy.sin(freestream.alpha - panels[-1].beta) )
    return b



def build_singularity_matrix(A_source, B_vortex):
    """
    Builds the left-hand side matrix of the system
    arising from source and vortex contributions.
    
    Parameters
    ----------
    A_source: 2D Numpy array of floats
        Source contribution matrix for the normal velocity.
    B_vortex: 2D Numpy array of floats
        Vortex contribution matrix for the normal velocity.
    
    Returns
    -------
    A:  2D Numpy array of floats
        Matrix of the linear system.
    """
    A = numpy.empty((A_source.shape[0] + 1, A_source.shape[1] + 1), dtype=float)
    # source contribution matrix
    A[:-1, :-1] = A_source
    # vortex contribution array
    A[:-1, -1] = numpy.sum(B_vortex, axis=1)
    # Kutta condition array
    A[-1, :] = kutta_condition(A_source, B_vortex)
    return A

def kutta_condition(A_source, B_vortex):
    """
    Builds the Kutta condition array.
    
    Parameters
    ----------
    A_source: 2D Numpy array of floats
        Source contribution matrix for the normal velocity.
    B_vortex: 2D Numpy array of floats
        Vortex contribution matrix for the normal velocity.
    
    Returns
    -------
    b: 1D Numpy array of floats
        The left-hand side of the Kutta-condition equation.
    """
    b = numpy.empty(A_source.shape[0] + 1, dtype=float)
    # matrix of source contribution on tangential velocity
    # is the same than
    # matrix of vortex contribution on normal velocity
    b[:-1] = B_vortex[0, :] + B_vortex[-1, :]
    # matrix of vortex contribution on tangential velocity
    # is the opposite of
    # matrix of source contribution on normal velocity
    b[-1] = - numpy.sum(A_source[0, :] + A_source[-1, :])
    return b

def vortex_contribution_normal(panels):
    """
    Builds the vortex contribution matrix for the normal velocity.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
    
    Returns
    -------
    A: 2D Numpy array of floats
        Vortex contribution matrix.
    """
    A = numpy.empty((panels.size, panels.size), dtype=float)
    # vortex contribution on a panel from itself
    numpy.fill_diagonal(A, 0.0)
    # vortex contribution on a panel from others
    for i, panel_i in enumerate(panels):
        for j, panel_j in enumerate(panels):
            if i != j:
                A[i, j] = -0.5 / numpy.pi * integral(panel_i.xc, panel_i.yc, 
                                                     panel_j,
                                                     numpy.sin(panel_i.beta),
                                                     -numpy.cos(panel_i.beta))
    return A

def source_contribution_normal(panels):
    """
    Builds the source contribution matrix for the normal velocity.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
    
    Returns
    -------
    A: 2D Numpy array of floats
        Source contribution matrix.
    """
    A = numpy.empty((panels.size, panels.size), dtype=float)
    # source contribution on a panel from itself
    numpy.fill_diagonal(A, 0.5)
    # source contribution on a panel from others
    for i, panel_i in enumerate(panels):
        for j, panel_j in enumerate(panels):
            if i != j:
                A[i, j] = 0.5 / numpy.pi * integral(panel_i.xc, panel_i.yc, 
                                                    panel_j,
                                                    numpy.cos(panel_i.beta),
                                                    numpy.sin(panel_i.beta))
    return A

def integral(x, y, panel, dxdk, dydk):
    """
    Evaluates the contribution from a panel at a given point.
    
    Parameters
    ----------
    x: float
        x-coordinate of the target point.
    y: float
        y-coordinate of the target point.
    panel: Panel object
        Panel whose contribution is evaluated.
    dxdk: float
        Value of the derivative of x in a certain direction.
    dydk: float
        Value of the derivative of y in a certain direction.
    
    Returns
    -------
    Contribution from the panel at a given point (x, y).
    """
    def integrand(s):
        return (((x - (panel.xa - numpy.sin(panel.beta) * s)) * dxdk +
                 (y - (panel.ya + numpy.cos(panel.beta) * s)) * dydk) /
                ((x - (panel.xa - numpy.sin(panel.beta) * s))**2 +
                 (y - (panel.ya + numpy.cos(panel.beta) * s))**2) )
    return integrate.quad(integrand, 0.0, panel.length)[0]

def vel_field(panels, freestream, X, Y):
    '''
    Computes the velocity field on a given 2D mesh.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
    freestream: Freestream object
        Freestream conditions.
    X: 2D Numpy array of floats
        x-coordinates of the mesh points.
    Y: 2D Numpy array of floats
        y-coordinates of the mesh points.
        
    Returns
    -------
    u: 2D Numpy array of floats
        x-component of the velocity vector field.
    v: 2D Numpy array of floats
        y-component of the velocity vector field.
        
        '''
    

    # freestream contribution
    u = freestream.u_inf * np.cos(freestream.alpha) * np.ones_like(X, dtype=float)
    v = freestream.u_inf * np.sin(freestream.alpha) * np.ones_like(X, dtype=float)

    # add the contribution from each source
    vec_integral = np.vectorize(integral)
    for panel in panels:
        u += panel.sigma / (2.0 * math.pi) * vec_integral(X, Y, panel, 1.0, 0.0)
        v += panel.sigma / (2.0 * math.pi) * vec_integral(X, Y, panel, 0.0, 1.0)

    return u, v


def define_panels(x, y, N=40, dir=1):
    """
    Discretizes the geometry into panels using 'cosine' method.
    
    Parameters
    ----------
    x: 1D array of floats
        x-coordinate of the points defining the geometry.
    y: 1D array of floats
        y-coordinate of the points defining the geometry.
    N: integer, optional
        Number of panels;
        default: 40.
    
    Returns
    -------
    panels: 1D Numpy array of Panel objects.
        The list of panels.
    """

    # If the direction is clockwise, reverse the order of the points
    # if dir<=0:
    #     x = x[::-1]
    #     y = y[::-1]
    
    R = (x.max() - x.min()) / 2.0  # circle radius
    x_center = (x.max() + x.min()) / 2.0  # x-coordinate of circle center
    
    theta = numpy.linspace(0.0, 2.0 * numpy.pi, N + 1)  # array of angles
    x_circle = x_center + R * numpy.cos(theta)  # x-coordinates of circle
    
    x_ends = numpy.copy(x_circle)  # x-coordinate of panels end-points
    y_ends = numpy.empty_like(x_ends)  # y-coordinate of panels end-points
    
    # extend coordinates to consider closed surface
    x, y = numpy.append(x, x[0]), numpy.append(y, y[0])
    
    # compute y-coordinate of end-points by projection
    I = 0
    for i in range(N):
        while I < len(x) - 1:
            if (x[I] <= x_ends[i] <= x[I + 1]) or (x[I + 1] <= x_ends[i] <= x[I]):
                break
            else:
                I += 1
        a = (y[I + 1] - y[I]) / (x[I + 1] - x[I])
        b = y[I + 1] - a * x[I + 1]
        y_ends[i] = a * x_ends[i] + b
    y_ends[N] = y_ends[0]
    
    # create panels
    panels = numpy.empty(N, dtype=object)
    for i in range(N):
        panels[i] = Panel(x_ends[i], y_ends[i], x_ends[i + 1], y_ends[i + 1])
    
    return panels




def main_pontoon_calc(foil_width, turbine_width, turbine_length, river_vel, plot=False):
    '''
    Complete function combining the potential flow solver to calculate the mean velocity between the pontoons
    Author: Phil Blecher. Built on framework by Lorena Barba. Barba, Lorena A., and Mesnard, Olivier (2019). Aero Python: classical aerodynamics of potential flow using Python. 
    Journal of Open Source Education, 2(15), 45, https://doi.org/10.21105/jose.00045

    Note: pontoon shape default is NACA 00xx

    Parameters
    ----------
    foil_width : float
        Width of the foil (m)
    turbine_width : float
        Width of the turbine (m) (between the inside edges of the pontoons)
    turbine_length : float
        Length of the turbine (m)
    river_vel : float
        Freestream river velocity (m/s)

    Returns
    -------
    vel : float
        Mean velocity between the pontoons (m/s)

    
    '''

    # create 2 pontoons
    x,y = naca_foil(foil_width)
    x2,y2 = naca_foil(foil_width)

    # scale the aerofoil length
    x = x * turbine_length
    x2 = x2 * turbine_length

    # generate the panels
    aerofoils = np.empty(2, dtype=object)
    coords = [[x,y],[x2,y2]]
    # discretise the aerofoil into panels
    for i,aerofoil in enumerate(coords):
        aerofoils[i] = define_panels(aerofoil[0],aerofoil[1], N=40, dir=i) # flip the direction in which the panels generated!

    # now move the aerofoils to the correct position
    for panel in aerofoils[1]:
        panel.ya = panel.ya + turbine_width + foil_width
        panel.yb = panel.yb + turbine_width + foil_width
        panel.yc = panel.yc + turbine_width + foil_width

    y2 = y2 + turbine_width + foil_width 

    # define the freestream conditions
    freestream = Freestream(river_vel,0) # freestream velocity, angle of attack

    # compute the source and vortex contribution matrices
    for i,p in enumerate(aerofoils):  
        A_source = source_contribution_normal(p)
        B_vortex = vortex_contribution_normal(p)

        A = build_singularity_matrix(A_source, B_vortex)
        b = build_freestream_rhs(p, freestream)

        # solve for singularity strengths
        strengths = numpy.linalg.solve(A, b)

        # store source strength on each panel
        for i , panel in enumerate(p):
            panel.sigma = strengths[i]
            
        # store circulation density
        gamma = strengths[-1]
        # tangential velocity at each panel center.
        compute_tangential_velocity(p, freestream, gamma, A_source, B_vortex)
        # surface pressure coefficient
        compute_pressure_coefficient(p, freestream)

    # compute the velocity field on the mesh grid
    # define velocity field
    nx, ny = 30, 30
    x_start, x_end = -0.3, turbine_length + 0.3
    y_start, y_end = - 0.3 - foil_width, turbine_width + 2*foil_width + 0.3
    x_ = np.linspace(x_start, x_end, nx)
    y_ = np.linspace(y_start, y_end, ny)
    X, Y = np.meshgrid(x_, y_)

    # compute the velocity field on the mesh grid
    u_tot = np.zeros((ny, nx))
    v_tot = np.zeros((ny, nx))

    u1, v1 = vel_field(aerofoils[0], freestream, X, Y)
    u2, v2 = vel_field(aerofoils[1], freestream, X, Y)

    u_tot = u2 + u1 
    v_tot = v2  + v1
    u_tot = u_tot / 2
    v_tot = v_tot / 2

    # compute the mean velocity between the pontoons
    # the area of interest is 0.3 of the way along the length of the turbine (the thickest part of the turbine)
    # and the mean of the velocities between [foil_width/2, turbine_width + foil_width/2]
    # in the net velocity field
    x_point = turbine_length*0.3
    y_points = [foil_width/2, turbine_width + foil_width/2]

    up_tot=0
    vp_tot=0
    num_points = 10 # number of points to average over
    for i,panels in enumerate(aerofoils):
        for y_point in np.linspace(y_points[0], y_points[1], num_points):
            u, v = vel_field(panels, freestream, x_point, y_point)
            up_tot += u
            vp_tot += v
    # normalize the velocity field
    up_tot /= 2*num_points 
    vp_tot /= 2*num_points

    # find the net velocity at this point
    vel = np.sqrt(up_tot**2 + vp_tot**2)
    

    if plot:
        plt.figure()
        plt.axis('equal')
        for panels in aerofoils:
            plt.fill([panel.xc for panel in panels],
                [panel.yc for panel in panels],
                color='k', linestyle='solid', linewidth=2, zorder=2)
            
        # add contours of velocity
        plt.contourf(X, Y, np.sqrt(u_tot**2 + v_tot**2), cmap='jet', levels=100)
        cbar = plt.colorbar(orientation='vertical', shrink=0.5, pad=0.1)
        cbar.set_label('Velocity magnitude', fontsize=16)

        plt.grid()
        plt.title('Flow around two pontoons', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.show()


    return vel