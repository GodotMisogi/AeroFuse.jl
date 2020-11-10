import numpy as np
from scipy import integrate
import math

class Source2D:
    """
    Defines a source or sink.
    
    Parameters
    ---
    sigma: strength
    x_0: source x-location
    y_0: source y-location
    """
    
    def __init__(self, sigma, x_0, y_0):
        self.strength = sigma
        self.x_0 = x_0
        self.y_0 = y_0
        
    def velocity(self, x, y):
        u = self.strength/(2*np.pi)*(x - self.x_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        v = self.strength/(2*np.pi)*(y - self.y_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        
        return (u, v)
    
    def stream_func(self, x, y):
        psi = self.strength/(2*np.pi)*np.arctan2(y - self.y_0, x - self.x_0)
        
        return psi
        
    def potential(self, x, y):
        pot = self.strength/(4*np.pi)*np.log(np.sqrt((x - self.x_0)**2 + (y - self.y_0)**2))
        
        return pot
    
class Uniform2D:
    """
    Defines a uniform flow.
    
    Parameters:
    1. u_inf = u
    2. alpha = angle
    """
    
    def __init__(self, u_inf=1.0, alpha=0.0):
        self.U = u_inf
        self.angle = np.radians(alpha)
        
    def velocity(self, x, y):
        u = self.U*(np.cos(self.angle))
        v = self.U*(np.sin(self.angle))
        
        return (u, v)
    
    def stream_func(self, x, y):
        psi = self.U*(y*np.cos(self.angle) - x*np.sin(self.angle))
        
        return psi
        
    def potential(self, x, y):
        pot = self.U*(x*np.cos(self.angle) + y*np.sin(self.angle))
        
        return pot

class Doublet2D:
    """
    Defines a doublet.
    
    Parameters
    ---
    kappa: strength
    x_0: doublet x-location
    y_0: doublet y-location
    """
    
    def __init__(self, kappa, x_0, y_0):
        self.strength = kappa
        self.x_0 = x_0
        self.y_0 = y_0
        
    def velocity(self, x, y):
        u = -self.strength/(2*np.pi)*((x - self.x_0)**2 - (y - self.y_0)**2)/((x - self.x_0)**2 + (y - self.y_0)**2)**2
        v = -self.strength/(2*np.pi)*2*(x - self.x_0)*(y - self.y_0)/((x - self.x_0)**2 + (y - self.y_0)**2)**2
        
        return (u, v)
    
    def stream_func(self, x, y):
        psi = -self.strength/(2*np.pi)*(y - self.y_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        
        return psi
        
    def potential(self, x, y):
        # CHECK THIS YOU RETARD
        pot = -self.strength/(2*np.pi)*(x - self.x_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        
        return pot

class Vortex2D:
    """
    Defines a vortex.
    
    Parameters
    ---
    Gamma: strength
    x_0: vortex x-location
    y_0: vortex y-location
    """
    
    def __init__(self, gamma, x_0, y_0):
        self.strength = gamma
        self.x_0 = x_0
        self.y_0 = y_0
        
    def velocity(self, x, y):
        u = -self.strength/(2*np.pi)*(y - self.y_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        v = self.strength/(2*np.pi)*(x - self.x_0)/((x - self.x_0)**2 + (y - self.y_0)**2)
        
        return (u, v)
    
    def stream_func(self, x, y):
        psi = -self.strength/(4*np.pi)*np.log(np.sqrt((x - self.x_0)**2 + (y - self.y_0)**2))
        
        return psi
        
    def potential(self, x, y):
        pot = self.strength/(2*np.pi)*np.arctan2(y - self.y_0, x - self.x_0)
        
        return pot

class InfiniteVortices:
    """
    Defines an infinite row of vortices.
    
    Parameters
    ---
    Gamma: strength
    a: spacing
    """
    
    def __init__(self, gamma, a):
        self.strength = gamma
        self.a = a
        
    def velocity(self, x, y):
        u = -self.strength/(2*self.a)*np.sinh(2*np.pi*y/self.a)/(np.cosh(2*np.pi*y/self.a) - np.cos(2*np.pi*x/self.a))
        v = self.strength/(2*self.a)*np.sin(2*np.pi*x/self.a)/(np.cosh(2*np.pi*y/self.a) - np.cos(2*np.pi*x/self.a))
        
        return (u, v)

def joukowski(z, c):
    """
    Computes the Joukowski transformation of a complex number z.
    
    Parameters:
    1. z: complex number in Cartesian coordinates
    2. c: some parameter
    """
    return z + c**2/z

def pressure_coefficient2D(u, v, u_inf):
    "Computes the pressure coefficient in 2 dimensions."
    cp = 1. - (u**2 + v**2)/u_inf**2
    return cp

class Panel2D:
    """
    Defines a panel in a given line segment.
    
    Parameters
    ---
    xs, ys: starting coordinates of the panel
    xe, ye: ending coordinates of the panel
    """
    
    def __init__(self, xs, ys, xe, ye):
        self.xs, self.ys = xs, ys
        self.xe, self.ye = xe, ye

        self.xc, self.yc = (xe + xs)/2.0, (ye + ys)/2.0
        self.length = np.sqrt((xe - xs)**2 + (ye - ys)**2)

        if xe - xs <= 0:
            self.angle = np.arccos((ye - ys)/self.length)
        else:
            self.angle = np.pi + np.arccos(-(ye - ys)/self.length)

        if self.angle <= np.pi:
            self.loc = 'upper'
        else:
            self.loc = 'lower'

class SourcePanel2D(Panel2D):
    """
    Modifies a panel by adding an infinite number of sources to it.
    """    

    def __init__(self, xs, ys, xe, ye):
        Panel2D.__init__(self, xs, ys, xe, ye)
        self.strength = 0.0 # Source strength
        self.vt = 0.0 # Tangential velocity
        self.cp = 0.0 # Pressure coefficient

    def potential(self, x, y):
        "Evaluates the potential at a point (x,y) for the panel."
        def integrand(s, x, y):
            return (np.log((x - (self.xs - np.sin(self.angle)*s))**2 + (y - (self.ys + np.cos(self.angle)*s))**2))
        
        pot = np.vectorize(lambda x, y: self.strength/(4*np.pi)*integrate.quad(integrand, 0.0, self.length, args=(x,y))[0])

        return pot(x,y)

    def integral(self, x, y, dxdz, dydz):
        "Evaluates the contribution of panel at a given location (x,y) of the current panel along a given direction z."
        def integrand(s):
            return (((x - (self.xs - np.sin(self.angle)*s))*dxdz 
                   + (y - (self.ys + np.cos(self.angle)*s))*dydz)
                   /((x - (self.xs - np.sin(self.angle)*s))**2 
                   + (y - (self.ys + np.cos(self.angle)*s))**2))

        return integrate.quad(integrand, 0.0, self.length)[0]

    def velocity(self, x, y):
        "Computes the velocity at a point (x,y) developed by the panel."        
        intregral = np.vectorize(self.integral)
        u = self.strength/(2*np.pi)*intregral(x, y, 1.0, 0.0)
        v = self.strength/(2*np.pi)*intregral(x, y, 0.0, 1.0)
        
        return u, v

def mag(xs):
    return np.sqrt(sum(map(lambda x: x**2, xs)))

# Transforms (x, y) to the coordinate system with (x_s, y_s) as origin oriented at angle_s.
def panelCoords(x, y, x_s, y_s, angle_s): 
    return rotation(x - x_s, y - y_s, angle_s)

# Rotation matrices
def invRotation(x, y, angle): 
    return (x*np.cos(angle) - y*np.sin(angle), x*np.sin(angle) + y*np.cos(angle))
def rotation(x, y, angle): 
    return (x*np.cos(angle) + y*np.sin(angle), -x*np.sin(angle) + y*np.cos(angle))

def pressureCoefficient2D(u, v, freestream_speed):
    return 1. - (u**2 + v**2)/freestream_speed**2

class DoubletPanel2D:

    def __init__(self, xs, ys, xe, ye, strength=0.0, vt=0.0, cp=0.0):
        self.xs, self.ys, self.xe, self.ye = xs, ys, xe, ye
        self.strength, self.vt, self.cp =  0, 0, 0
        self.length = mag([ xe - xs, ye - ys ])
        self.angle = np.pi - np.arctan2(ye - ys, xe - xs)
        self.loc = "lower" if (np.pi <= self.angle <= 2*np.pi) else "upper"
        self.normal = np.array([np.cos(self.angle), np.sin(self.angle)])
        self.tangent = np.array([-np.sin(self.angle), np.cos(self.angle)])
    
        # Collocation points in global coordinates
        self.xc, self.yc = (xe + xs)/2.0, (ye + ys)/2.0

        # Panel coordinates
        self.xsl, self.ysl, self.xel, self.yel  = 0., 0., self.length, 0.
        self.xcl, self.ycl = self.length/2., 0

    def influence(self, xg, yg):
        x, y = panelCoords(xg, yg, self.xs, self.ys, self.angle)
        return [1/(2*np.pi)*(y/((x - self.xsl)**2 + y**2) - y/((x - self.xel)**2 + y**2)), -1/(2*np.pi)*((x - self.xsl)/((x - self.xsl)**2 + y**2) - (x - self.xel)/((x - self.xel)**2 + y**2))]      

    def potential(self, x, y):
        return self.strength/(4*np.pi)*(x - self.xc)/((x - self.xc)**2 + (y - self.yc)**2)

    def velocity(self, xg, yg): 
        x, y = panelCoords(xg, yg, self.xs, self.ys, self.angle)
        return [self.strength/(2*np.pi)*(y/((x - self.xsl)**2 + y**2) - y/((x - self.xel)**2 + y**2)), -self.strength/(2*np.pi)*((x - self.xsl)/((x - self.xsl)**2 + y**2) - (x - self.xel)/((x - self.xel)**2 + y**2))]

class DoubletPanelSolver2D:
    """
    Applies the doublet panel method to a list of panels and a uniform flow.
    """
    def __init__(self, panels, uniform):
        self.panels = panels
        self.num_panels = len(panels)
        self.woke_panel = DoubletPanel2D(self.panels[-1].xe, self.panels[-1].ye, 1000*self.panels[-1].xe, self.panels[-1].ye)
        self.uniform = uniform
        self.u = [uniform.U*np.cos(uniform.angle), uniform.U*np.sin(uniform.angle)]
        
    def solveStrengths(self):
        """
        Solves for the doublet strengths of all panels.
        """
        # Construct doublet influence matrix for normal direction.
        doublet = np.array([ [ -2/(np.pi*panel_i.length) if i == j else np.sum(np.multiply(panel_j.influence(panel_i.xc, panel_i.yc), rotation(*panel_i.normal, panel_j.angle))) for (i, panel_i) in enumerate(self.panels) ] for (j, panel_j) in enumerate(self.panels) ])

        woke = [ sum(np.multiply(self.woke_panel.influence(panel.xc, panel.yc), rotation(*panel.normal, self.woke_panel.angle))) for panel in self.panels ]

        # Kutta condition
        kutta = np.zeros(self.num_panels + 1)
        kutta[0] = 1.
        kutta[-2] = -1.
        kutta[-1] = 1.

        # Matrix construction
        An = np.zeros((self.num_panels + 1, self.num_panels + 1))
        An[:-1, :-1] = doublet
        An[:-1, -1] = woke
        An[-1, :] = kutta

        # Construct freestream RHS.
        bn = np.append([ np.sum(np.multiply(self.u, panel.normal)) for panel in self.panels ], 0)
        
        # Solve system.
        strengths = np.linalg.solve(An, bn)

        # Update panel strengths.
        for (panel, strength) in zip(self.panels, strengths):
            panel.strength = strength
        return strengths
    
    def aerodynamicsss(self, pressure = True):
        """
        Solves for the velocities of all panels and their pressure coefficients.
        """
        pans = np.append(self.panels, self.woke_panel)

        # Construct matrix for tangential direction.
        At = np.array([ [ 0. if i == j else np.sum(np.append(panel_i.velocity(panel_j.xc, panel_j.yc), panel_j.tangent)) for (i, panel_i) in enumerate(pans)] for (j, panel_j) in enumerate(pans) ])

        # Surface velocities
        qts = np.sum(At, axis=1)
        # print(qts)

        # Construct freestream tangential condition
        bt = [ np.sum(np.append(self.u, panel.tangent)) for panel in pans ]

        # Solve system.
        vts = qts + bt

        # Update panel velocities and pressure coefficients.
        for (panel, vt) in zip(pans, vts):
            panel.vt = vt
            panel.cp = pressureCoefficient2D(0., vt, self.uniform.U) if pressure else 0
        return vts

    def liftCoefficient(self):
        xs = [ panel.xs for panel in self.panels ]
        c = abs(max(xs) - min(xs))
        pans = np.append(self.panels, self.woke_panel)
        strengths = [ panel.strength for panel in pans ]
        strengths1 = [ strengths[i+1] - strengths[i] for i in range(len(strengths) - 1) ]
        cl = (-2/self.uniform.U*np.sum(strengths1), 
        -2*strengths[-1]/self.uniform.U, 
        2/self.uniform.U*np.sum([ -panel1.cp*mag([panel1.xc - panel2.xc, panel1.yc - panel2.yc])*np.cos(panel1.angle)/c for (panel1, panel2) in zip(pans[:-1], pans[1:])]))

        return cl
    
    # Compute errors
    def error(self):
        error = sum([ panel.strength*panel.length for panel in self.panels ])
        print(f'Sum of sources: {error}')

        return error

    # Compute velocities
    def velocity(self, x, y):
        "Computes the velocity at a point (x, y) generated by all the panels."
        vels = [ self.uniform.velocity(x, y)[i] + sum([ panel.velocity(x, y)[i] for panel in self.panels ]) for i in range(2) ]
         
        return vel

    def potential(self, x, y):
        "Computes the potential at a point (x, y) generated by all the panels."
        return self.uniform.potential(x, y) + sum([ panel.potential(x, y) for panel in self.panels ])

        

class SourcePanel2DSolver():
    """
    Applies the source panel method to a list of panels and a uniform flow.
    """

    def __init__(self, panels, uniform2D):
        self.panels = panels
        self.freestream = uniform2D

    def solve_sources(self):
        """
        Solves for the source strengths of all panels.
        """
        # Construct matrix for normal direction.
        An = np.array([ [ 0.5 if i == j else 0.5/np.pi*p_i.integral(p_j.xc, p_j.yc, np.cos(p_j.angle), np.sin(p_j.angle)) for i, p_i in enumerate(self.panels) ] for j, p_j in enumerate(self.panels) ], dtype=float)
        
        # Construct freestream RHS.
        bn = -self.freestream.U*np.cos([ self.freestream.angle - panel.angle for panel in self.panels ]) 

        # Solve system.
        sources = np.linalg.solve(An, bn)

        # Update panel strengths.
        for panel, strength in zip(self.panels, sources):
            panel.strength = strength
        
        return sources

    def error(self):
        error = sum([ panel.strength*panel.length for panel in self.panels ])
        print(f'Sum of sources: {error}')

        return error

    def tangential_velocity(self, pressure=True):
        """
        Solves for the velocities of all panels and their pressure coefficients.
        """
        # Construct matrix for tangential direction.
        At = np.array([ [ 0.0 if i == j else 0.5/np.pi*p_i.integral(p_j.xc, p_j.yc, -np.sin(p_j.angle), np.cos(p_j.angle)) for i, p_i in enumerate(self.panels) ] for j, p_j in enumerate(self.panels) ], dtype=float)

        # Construct freestream RHS.
        bt = self.freestream.U*np.sin([ self.freestream.angle - panel.angle for panel in self.panels ]) 

        # Solve system.
        vts = np.dot(At, [ panel.strength for panel in self.panels ] ) + bt

        # Update panel velocities and pressure coefficients.
        for panel, vt in zip(self.panels, vts):
            panel.vt = vt
            if pressure:
                panel.cp = pressure_coefficient2D(0, vt, self.freestream.U)

        return vts

    def velocity(self, x, y):
        "Computes the velocity at a point (x, y) generated by all the panels."
        vels = [ self.freestream.velocity(x, y)[i] + sum([ panel.velocity(x, y)[i] for panel in self.panels ]) for i in range(2) ]
         
        return vels

    def potential(self, x, y):
        "Computes the potential at a point (x, y) generated by all the panels."
        return self.freestream.potential(x, y) + sum([ panel.potential(x, y) for panel in self.panels ])

class VortexSourcePanel2DSolver():
    """
    Applies the vortex source panel method to a list of panels and a uniform flow.
    """

    def __init__(self, panels, uniform2D):
        self.panels = panels
        self.freestream = uniform2D
        self.gamma = 0.0 # Vortex strength for panels
        
        # Construct source matrix for normal direction.
        self.source_matrix = np.array([ [ 0.5 if i == j else 0.5/np.pi*p_i.integral(p_j.xc, p_j.yc, np.cos(p_j.angle), np.sin(p_j.angle)) for i, p_i in enumerate(self.panels) ] for j, p_j in enumerate(self.panels) ], dtype=float)

        # Vortex matrix
        self.vortex_matrix = np.array([ [ 0.0 if i == j else -0.5/np.pi*p_i.integral(p_j.xc, p_j.yc, np.sin(p_j.angle), -np.cos(p_j.angle)) for i, p_i in enumerate(self.panels) ] for j, p_j in enumerate(self.panels) ], dtype=float)

        self.num_panels = len(self.panels)

    def solve_strengths(self):
        """
        Solves for the source strengths of all panels.
        """
        
        # Kutta condition
        kutta = np.empty(self.num_panels + 1, dtype=float)
        kutta[:-1] = self.vortex_matrix[0, :] + self.vortex_matrix[-1, :]
        kutta[-1] = -sum(self.source_matrix[0, :] + self.source_matrix[-1, :]) 

        # Global matrix construction
        An = np.empty((self.num_panels + 1, self.num_panels + 1), dtype=float)
        An[:-1, :-1] = self.source_matrix
        An[:-1, -1] = [ sum(row) for row in self.vortex_matrix ]
        An[-1, :] = kutta

        # Construct freestream RHS.
        bn = np.append(-self.freestream.U*np.cos([self.freestream.angle - panel.angle for panel in self.panels ]), -self.freestream.U * (np.sin(self.freestream.angle - self.panels[0].angle) + np.sin(self.freestream.angle - self.panels[-1].angle)))

        # Solve system.
        strengths = np.linalg.solve(An, bn)

        # Update panel strengths.
        for panel, strength in zip(self.panels, strengths):
            panel.strength = strength

        # Update vorticity.
        self.gamma = strengths[-1]

        return strengths

    def error(self):
        error = sum([ panel.strength*panel.length for panel in self.panels ])
        print(f'Sum of sources: {error}')

        return error

    def tangential_velocity(self, pressure=True):
        """
        Solves for the velocities of all panels and their pressure coefficients.
        """
        At = np.empty((self.num_panels, self.num_panels + 1), dtype=float)

        # Construct matrix for tangential direction.
        At[:, :-1] = self.vortex_matrix
        At[:, -1] = [ -sum(row) for row in self.source_matrix ]
        
        # Construct freestream RHS.
        bt = self.freestream.U*np.sin([ self.freestream.angle - panel.angle for panel in self.panels ]) 

        # Solve system.
        vts = np.dot(At, np.append([ panel.strength for panel in self.panels], self.gamma) ) + bt

        # Update panel velocities and pressure coefficients.
        for panel, vt in zip(self.panels, vts):
            panel.vt = vt
            if pressure:
                panel.cp = pressure_coefficient2D(0, vt, self.freestream.U)

        return vts

    def velocity(self, x, y):
        "Computes the velocity at a point (x, y) generated by all the panels."
        vels = [ self.freestream.velocity(x, y)[i] + sum([ panel.velocity(x, y)[i] for panel in self.panels ]) for i in range(2) ]
         
        return vels

    def potential(self, x, y):
        "Computes the potential at a point (x, y) generated by all the panels."
        return self.freestream.potential(x, y) + sum([ panel.potential(x, y) for panel in self.panels ])

    def lift_coefficient(self):
        # Lift coefficient
        c = abs(max(panel.xs for panel in self.panels) -
                min(panel.xs for panel in self.panels))
        cl = (self.gamma * sum(panel.length for panel in self.panels) /
            (0.5 * self.freestream.U * c))
        print(f'Lift coefficient: CL = {cl}')

        return cl

def cosine_panels(x, y, n = 40):
    """
    Discretises a geometry consisting of x and y coordinates into panels by projecting the x-coordinate of a circle onto the geometry.
    """
    r = (x.max() - x.min())/2
    x_center = (x.max() + x.min())/2
    x_circ = x_center + r*np.cos(np.linspace(0.0, 2*np.pi, n + 1))

    x_ends = np.copy(x_circ)
    y_ends = np.empty_like(x_ends)

    x, y = np.append(x, x[0]), np.append(y, y[0])

    j = 0
    for i in range(n):
        while j < len(x) - 1:
            if (x[j] <= x_ends[i] <= x[j+1]) or (x[j+1] <= x_ends[i] <= x[j]):
                break
            else:
                j += 1
        m = (y[j+1] - y[j])/(x[j+1] - x[j])
        c = y[j+1] - m*x[j+1]
        y_ends[i] = m*x_ends[i] + c
    y_ends[n] = y_ends[0]
    
    panels = [ Panel2D(x_ends[i], y_ends[i], x_ends[i+1], y_ends[i+1]) for i in range(n) ]

    return panels
        