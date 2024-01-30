import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import integrate

'''
Simulate the flow around two aerofoils side by side, to simulate a cross-section of a pontoon

The aerofoils are identical and in water, so the flow is incompressible and inviscid

The flow is assumed to be 2D, so the flow is in the x-y plane

Start point: Barba, Lorena A., and Mesnard, Olivier (2019). Aero Python: classical aerodynamics of potential flow using Python. 
Journal of Open Source Education, 2(15), 45, https://doi.org/10.21105/jose.00045

Using source panel method
'''

'''
Aerofoil geometries:

'''
def nozzle_foil(alpha=0):
    '''
    a 1d line for each side of the nozzle
    
    '''

    x = np.linspace(0, 1, 100)
    # y is a diagonal line at angle alpha
    y = x*np.tan(np.radians(alpha))

    return x, y


def custom_foil(t, t_c ,alpha=0, side='left'):
    '''
    Left (LE) to right (TE) will:
    start with a sideways parabolic until reaching a width of t/2 at x=t_c,
    after which will be flat to the TE
    the other side will follow the NACA shape

    Parameters
    ----------
    t: float
        maximum thickness of the aerofoil as a fraction of the chord length
        (so NACA0012 has t = 0.12)
    t_c: float
        x value at which the maximum thickness occurs on the inside edge
    alpha: float
        angle of attack deg
    side: str
        left or right pontoon
    
    '''

    if side == 'left':
        # define the x coordinates of the upper surface
        x_upper_1 = np.linspace(0, t_c, 50)
        x_upper_2 = np.linspace(t_c,1,50)
        # define the y coordinates of the upper surface
        y_upper_1 = 5*t * (0.2969 * np.sqrt(x_upper_1) - 0.1260 * x_upper_1 - 0.3516 * x_upper_1**2 + 0.2843 * x_upper_1**3 - 0.1015 * x_upper_1**4)
        # for the second part will be diagonal down to the TE:
        y_upper_2 = np.ones(len(x_upper_2))*y_upper_1[-1] - (y_upper_1[-1]/(1-t_c))*(x_upper_2-t_c)
        # combine the two parts
        x_upper = np.concatenate((x_upper_1, x_upper_2))
        y_upper = np.concatenate((y_upper_1, y_upper_2))
        # ref NACA
        # define the x coordinates of the lower surface
        x_lower = np.linspace(0, 1, 100)
        # define the y coordinates of the lower surface
        y_lower = -5*t * (0.2969 * np.sqrt(x_lower) - 0.1260 * x_lower - 0.3516 * x_lower**2 + 0.2843 * x_lower**3 - 0.1015 * x_lower**4)
    elif side == 'right':
        # define the x coordinates of the upper surface
        x_lower_1 = np.linspace(0, t_c, 50)
        x_lower_2 = np.linspace(t_c,1,50)
        # define the y coordinates of the upper surface
        y_lower_1 = -5*t * (0.2969 * np.sqrt(x_lower_1) - 0.1260 * x_lower_1 - 0.3516 * x_lower_1**2 + 0.2843 * x_lower_1**3 - 0.1015 * x_lower_1**4)
        # for the second part will be flat diagonal to the TE:
        y_lower_2 = np.ones(len(x_lower_2))*y_lower_1[-1] - (y_lower_1[-1]/(1-t_c))*(x_lower_2-t_c)
        # combine the two parts
        x_lower = np.concatenate((x_lower_1, x_lower_2))
        y_lower = np.concatenate((y_lower_1, y_lower_2))
        # ref NACA
        # define the x coordinates of the lower surface
        x_upper = np.linspace(0, 1, 100)
        # define the y coordinates of the lower surface
        y_upper = 5*t * (0.2969 * np.sqrt(x_upper) - 0.1260 * x_upper - 0.3516 * x_upper**2 + 0.2843 * x_upper**3 - 0.1015 * x_upper**4)

    else:
        raise ValueError('define the side')

    # make it all one surface anticlockwise
    x = np.concatenate((x_upper[::-1], x_lower))
    y = np.concatenate((y_upper[::-1], y_lower))

    x_rot, y_rot = rotate(x,y,alpha=alpha)

    return x_rot, y_rot





def naca_one_side(t, alpha=0, side='left'):
    '''
    function returns the x and y coordinates of the aerofoil geometry

    Parameters
    ----------
    t: float
        maximum thickness of the aerofoil as a fraction of the chord length
        (so NACA0012 has t = 0.12)
    alpha: float
        angle of attack deg
    side: str
        left or right pontoon

    '''
    if side == 'left':
        # define the x coordinates of the upper surface
        x_upper = np.linspace(0, 1, 100)
        # define the y coordinates of the upper surface
        y_upper = 5*t * (0.2969 * np.sqrt(x_upper) - 0.1260 * x_upper - 0.3516 * x_upper**2 + 0.2843 * x_upper**3 - 0.1015 * x_upper**4)
        # ref NACA
        # define the x coordinates of the lower surface
        x_lower = np.linspace(0, 1, 100)
        # define the y coordinates of the lower surface
        y_lower = np.zeros(len(x_lower))
    elif side == 'right':
        # define the x coordinates of the upper surface
        x_lower = np.linspace(0, 1, 100)
        # define the y coordinates of the lower surface
        y_lower = -5*t * (0.2969 * np.sqrt(x_lower) - 0.1260 * x_lower - 0.3516 * x_lower**2 + 0.2843 * x_lower**3 - 0.1015 * x_lower**4)
        # ref NACA
        # define the x coordinates of the lower surface
        x_upper = np.linspace(0, 1, 100)
        # define the y coordinates of the lower surface
        y_upper = np.zeros(len(x_upper))

    else:
        raise ValueError('define the side')

    # make it all one surface anticlockwise
    x = np.concatenate((x_upper[::-1], x_lower))
    y = np.concatenate((y_upper[::-1], y_lower))

    x_rot, y_rot = rotate(x,y,alpha=alpha)

    return x_rot, y_rot


# Define the geometry of the aerofoil for symmetric aerofoil
def naca_foil(t, alpha=0):
    '''
    function returns the x and y coordinates of the aerofoil geometry

    Parameters
    ----------
    t: float
        maximum thickness of the aerofoil as a fraction of the chord length
        (so NACA0012 has t = 0.12)
    alpha: deg float
        angle of attack
    
    '''

    # define the x coordinates of the upper surface
    x_upper = np.linspace(0, 1, 100)
    # define the y coordinates of the upper surface
    y_upper = 5*t * (0.2969 * np.sqrt(x_upper) - 0.1260 * x_upper - 0.3516 * x_upper**2 + 0.2843 * x_upper**3 - 0.1015 * x_upper**4)
    # ref NACA
    # define the x coordinates of the lower surface
    x_lower = np.linspace(0, 1, 100)
    # define the y coordinates of the lower surface
    y_lower = -5*t * (0.2969 * np.sqrt(x_lower) - 0.1260 * x_lower - 0.3516 * x_lower**2 + 0.2843 * x_lower**3 - 0.1015 * x_lower**4)
    
    # make it all one surface anticlockwise
    x = np.concatenate((x_upper[::-1], x_lower))
    y = np.concatenate((y_upper[::-1], y_lower))

    x_rot, y_rot = rotate(x,y,alpha=alpha)

    return x_rot, y_rot

############################################################################################################
'''

Calculation functions and classes

'''

def rotate(x, y, alpha=0):

    if alpha>0:
        alpha = alpha-360

    alpha = np.radians(alpha)
    # rotate the aerofoil by alpha
    x_rot = x*np.cos(alpha) - y*np.sin(alpha)
    y_rot = x*np.sin(alpha) + y*np.cos(alpha)

    return x_rot, y_rot

# Define the panels
def define_panels(x, y, N = 20):
    '''
    Discretizes the geometry into panels using the 'cosine' method.

    Parameters
    ----------
    x: 1D array of floats
        x-coordinate of the points defining the geometry
    y: 1D array of floats
        y-coordinate of the points defining the geometry
    N: integer, optional
        Number of panels;
        default: 20
    alpha: float
        angle of attack
    
    '''

    R = (x.max() - x.min()) / 2  # radius of the circle
    x_center = (x.max() + x.min()) / 2  # x-coord of the center
    # define x-coord of the circle points
    x_circle = x_center + R * np.cos(np.linspace(0, 2 * math.pi, N + 1))

    x_ends = np.copy(x_circle)  # projection of the x-coord on the surface
    y_ends = np.empty_like(x_ends)  # initialization of the y-coord np array

    x, y = np.append(x, x[0]), np.append(y, y[0])  # extend arrays using np.append

    # computes the y-coordinate of end-points
    I = 0
    for i in range(N):
        while I < len(x) - 2:
            if (x[I] <= x_ends[i] <= x[I + 1]) or (x[I + 1] <= x_ends[i] <= x[I]):
                break
            else:
                I += 1
        a = (y[I + 1] - y[I]) / (x[I + 1] - x[I])
        b = y[I + 1] - a * x[I + 1]
        y_ends[i] = a * x_ends[i] + b
    y_ends[N] = y_ends[0]
    
    panels = np.empty(N, dtype=object)

    # define each panel
    for i in range(N):
        panels[i] = Panel(x_ends[i], y_ends[i], x_ends[i + 1], y_ends[i + 1])
    
    return panels

def integral(x, y, panel, dxdz, dydz):
    '''
    Evaluates the contribution of a panel at one point.

    Parameters
    ----------
    x: float
        x-coordinate of the target point
    y: float
        y-coordinate of the target point
    panel: Panel object
        Panel whose contribution is evaluated
    dxdz: float
        Derivative of x in the z-direction
    dydz: float
        Derivative of y in the z-direction
    
    Returns
    -------
    Integral over the panel of the influence at the target point
    
    '''
    def integrand(s):
        return (((x - (panel.xa - np.sin(panel.beta) * s)) * dxdz +
                (y - (panel.ya + np.cos(panel.beta) * s)) * dydz) /
                ((x - (panel.xa - np.sin(panel.beta) * s))**2 +
                (y - (panel.ya + np.cos(panel.beta) * s))**2))
    
    return integrate.quad(integrand, 0.0, panel.length)[0]


def A_mat(panels):
    '''
    Builds the source matrix.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
        
    Returns
    -------
    A: 2D Numpy array of floats
        Source matrix.
            
    '''
    N = len(panels)
    A = np.empty((N, N), dtype=float)
    # source contribution on a panel from itself
    np.fill_diagonal(A, 0.5)
    # source contribution on a panel from others
    for i, p_i in enumerate(panels):
        for j, p_j in enumerate(panels):
            if i != j:
                A[i, j] = 0.5 / math.pi * integral(p_i.xc, p_i.yc, p_j,
                                                np.cos(p_i.beta),
                                                np.sin(p_i.beta))
                
    return A


def b_vec(panels, freestream):
    '''
    Builds the RHS of the linear system.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
    freestream: Freestream object
        Freestream conditions.
        
    Returns
    -------
    b: 1D Numpy array of floats
        RHS of the linear system.
            
    '''
    N = len(panels)
    b = np.empty(N, dtype=float)
    # freestream contribution on a panel
    for i, panel in enumerate(panels):
        b[i] = -freestream.u_inf * np.cos(freestream.alpha - panel.beta)
        
    return b

def tan_vel(panels, freestream):
    '''
    Computes the tangential velocity on the surface of the panels.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
    freestream: Freestream object
        Freestream conditions.
        
        '''

    N = len(panels)
    A = np.empty((N, N), dtype=float)
    # source contribution on a panel from itself
    np.fill_diagonal(A, 0.0)
    # source contribution on a panel from others
    for i, p_i in enumerate(panels):
        for j, p_j in enumerate(panels):
            if i != j:
                A[i, j] = -0.5 / math.pi * integral(p_i.xc, p_i.yc, p_j,
                                                np.sin(p_i.beta),
                                                -np.cos(p_i.beta))
                
    b = freestream.u_inf * np.sin([freestream.alpha - panel.beta for panel in panels])

    sigma = np.array([panel.sigma for panel in panels])
    vt = np.dot(A, sigma) + b

    for i, panel in enumerate(panels):
        panel.vt = vt[i]

    # apply the kutta condition at the trailing and leading edge
    panels[0].vt = 0.0
    panels[N // 2].vt = 0.0
    panels[-1].vt = 0.0

    # apply no-slip condition at the trailing and leading edge
    panels[0].sigma = 0.0
    panels[N // 2].sigma = 0.0
    panels[-1].sigma = 0.0


def cp(panels, freestream):
    '''
    Computes the pressure coefficient on the surface of the panels.
    
    Parameters
    ----------
    panels: 1D array of Panel objects
        List of panels.
    freestream: Freestream object
        Freestream conditions.
        
        '''

    for panel in panels:
        panel.cp = 1.0 - (panel.vt / freestream.u_inf)**2

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


class Freestream:

    def __init__(self, u_inf = 1.0, alpha = 0.0):
        '''
        Sets the freestream conditions.

        Parameters
        ----------
        u_inf: float, optional
            Freestream speed;
            default: 1.0
        alpha: float, optional
            Angle of attack in degrees;
            default: 0.0
        
        '''
        self.u_inf = u_inf
        self.alpha = np.radians(alpha)  # degrees to radians


class Panel:

    def __init__(self, xa, ya, xb, yb):
        '''
        Initializes the panel.

        Parameters
        ----------
        xa: float
            x-coordinate of the first end-point
        ya: float
            y-coordinate of the first end-point
        xb: float
            x-coordinate of the second end-point
        yb: float
            y-coordinate of the second end-point
        
        '''
        self.xa, self.ya = xa, ya  # panel starting point
        self.xb, self.yb = xb, yb  # panel ending point

        self.xc, self.yc = (xa + xb) / 2, (ya + yb) / 2  # panel centre
        self.length = np.sqrt((xb - xa) ** 2 + (yb - ya) ** 2)  # panel length

        # orientation of panel
        if xb - xa <= 0.0:
            self.beta = np.arccos((yb - ya) / self.length)
        elif xb - xa > 0.0:
            self.beta = np.pi + np.arccos(-(yb - ya) / self.length)

        # location of panel
        if self.beta <= np.pi:
            self.loc = 'upper'
        else:
            self.loc = 'lower'

        self.sigma = 0.0  # source strength
        self.vt = 0.0  # tangential velocity
        self.cp = 0.0  # pressure coefficient

    def rotate_panel(self, alpha=0.0):
        self.xa, self.ya = rotate(self.xa,self.ya, alpha=alpha)
        self.xb, self.yb = rotate(self.xb,self.yb, alpha=alpha)