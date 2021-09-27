import numpy as np 
from scipy.integrate import solve_ivp
from .bisection import bisection
from .constants import  *

class CompactStar (object):

    def __init__ (self, eos):
        """
        Given an eos return mass and radius of the star (given also an external central pressure).

        Parameters
        ==========
        eos: callable function
            polytropicEOS() instance
        radii = array of float
            range of radii where to search the real radius of the star via tov eqs.
            this range must contain the expected radius of the star: for a ns we must have 
            100 km (1e5 m) as maximum radius, while for a white dwarf we shold have 5e7 km.
            range expressed in geom units, i.e in meters. Default value set to work with Neutron Stars.
        """       
        self.radii = np.arange(1e-3,1e5,0.1) # in m (geom); default value
        self.eos = eos

    def set_radii_range (self, radii):
        """
        Parameters
        ==========
        radii = array of float
            range of radii where to integrate structure eqs. radii are in meters

        NOTES:  we expect neutron stars to be within 50 km of radius and white dwarf much larger
                first value must be closed to zero to well define initial conditions
        """ 
        self.radii = radii

    def tov_eqs(self, r, y):
        """
        Function that returns TOV structure eqs for a Neutron Star. 
        Use of geometrized units system: G=c=1
        These eqs are written in the form dy/dx = f(x,y).

        Parameters
        ==========
        r : float
            independent variable (distance from the centre of the star)
        y : array of float
            dependent variable which is an array containing the mass M = y[0] and the pressure P = y[1]
        ----------
        dy : array of float
            array containing the TOV derivatives of M and of P, i.e. f(x, y)
        ----------
        NOTES: eden : obtained from sels.eos.eden_from_pressure(y[1])
        """
        # some clever name
        m, p = y

        # Call the EOS
        eden = self.eos.eden_from_pressure(p)
        # Set the RHS
        dy = np.zeros(len(y))
        dy[0] = 4*pi*eden*r**2
        dy[1] = -(eden+p)*(m + 4*pi*r**3*p)/(r*(r-2*m))
        return dy

    def newton_eqs(self, r, y):
        """
        Function that returns Newtonian structure eqs for a Neutron Star. 
        Use of geometrized units system: G=c=1
        These eqs are written in the form dy/dx = f(x,y).

        Parameters
        ==========
        r : float
            independent variable (distance from the centre of the star)
        y : array of float
            dependent variable which is an array containing the mass M = y[0] and the pressure P = y[1]
        ----------
        dy : array of float
            array containing the TOV derivatives of M and of P, i.e. f(x, y)
        ----------
        NOTES: eden : obtained from sels.eos.eden_from_pressure(y[1])
        """

        # some clever name
        m, p = y

        # Call the EOS
        eden = self.eos.eden_from_pressure(p)

        # Set the RHS
        dy = np.zeros(len(y))
        dy[1] = - (eden*m)/(r**2)
        dy[0] = 4.0*pi*r**2*eden

        return dy
  
    def set_initial_conditions(self, central_pressure): 
        """
        Utility to set initial conditions.

        Parameters:
        ===========
        - central_pressure : float
            pressure at r = r_min
        -----------
        - m0, p0: float
                initial mass and central pressure
        """
        p0 = central_pressure
        e0 = self.eos.eden_from_pressure(p0)
        m0 = 4.*pi*e0*self.radii[0]**3 
        return m0, p0

    def structure_solver(self, eqs_type, central_pressure, max_step = None):
        """
        Solve structure equations given a central pressure.

        Parameters:
        ===========
        p0 : float
            central pressure
        eqs_type : string
            specify if we solve TOV or Newton eqs. Must be 'TOV' or 'Newton'
        -----------
        r_out, p_out, m_out : array of floats
            p(r) and m(r) of the star with the given central pressure. 
            Last values of these arrays are the Radius, the Mass and the pressure at the surface of the star
        """

        eqs_dict = { 'Newton' : self.newton_eqs,
                     'TOV' : self.tov_eqs
                }

        if (eqs_type not in eqs_dict):
            raise ValueError("Must choose 'Newton' or 'TOV' in 'eqs_type'")
            quit()

        # we define now an attribute of CompactStars that should be passed to found_radius 
        #(found_radius il called only inside this function so it's ok)
        
        if self.eos.type == 'ImplicitEos' : 
            # we change the default fermi energy momenta range for eos in order to surely include all the pressures
            pressure_shifted = lambda x : self.eos.pressure_fermi_momentum(x) - central_pressure
            x_p = bisection(pressure_shifted, 0, 1e3, int(1e3))
            x_new = np.linspace(-int(x_p)-1, int(x_p)+1, 1000)
            self.eos.set_x_range(x_new)

        initial = self.set_initial_conditions(central_pressure)
        r_span = [self.radii[0], self.radii[-1]]
        
        if max_step is not None:
            max_step = max_step
        else: max_step = np.inf # default value of solve_ivp

        # we use a lambda function trick to pass central_pressure to 'events'
        # we set the 'terminal' attribute of found_radius to be True to stop integration properly
       
        stop = lambda r,y : found_radius(r, y, central_pressure)
        stop.terminal = True

        solutions = solve_ivp(eqs_dict[eqs_type], r_span, initial, method = 'RK45', max_step = max_step,
                events = stop)

        r_out = solutions.t *conversion_dict['geom']['lenght']['km']
        m_out = solutions.y[0] *conversion_dict['geom']['mass']['m_sol'] 
        p_out = solutions.y[1] *conversion_dict['geom']['pressure']['cgs'] 
            
        return r_out, m_out, p_out

    def mass_vs_radius (self, eqs_type, central_pressures, max_step = None):
        """
        Returns arrays containing the Mass and the Radius of the compact stars, given a range of central pressures.
        A couple of values (R,M) for eache pressure

        Parameters:
        ===========
        central_pressures : array of float
            central pressures
        eqs_type : string
            specify if we solve TOV or Newton eqs. Must be 'TOV' or 'Newton'
        -----------
        R_out, M_out : array of floats
            (R,M) of the star for each given central pressure. 
        """
        #print("Solving ", eqs_type, " equations...")
        final_pressure = central_pressures[-1]
        
        if self.eos.type == 'ImplicitEos' : 
            # we change the default fermi energy momenta range for eos in order to surely include all the pressures
            pressure_shifted = lambda x : self.eos.pressure_fermi_momentum(x) - final_pressure
            x_p = bisection(pressure_shifted, 0, 1e3, int(1e4))
            x_new = np.linspace(-int(x_p)-1, int(x_p)+1, 10000)
            self.eos.set_x_range(x_new)

        R_star = np.empty_like(central_pressures)
        M_star = np.empty_like(central_pressures)

        for i, pc in enumerate(central_pressures):
            r_star , m_star, p_star = self.structure_solver(eqs_type, pc, max_step)
            R_star[i] = r_star[-1]
            M_star[i] = m_star[-1]

        return R_star, M_star

# END OF CompactStar class
#===================================================================================#

# External function used to find when pressure goes under a certain value

def found_radius(t, y, central_pressure):
    """
    Search when the pressure goes under a certain default value wrt the central pressure.

    Parameters:
    ===========
    t : float
        independent variable. 
        Not used, but necessary to have the right signature required by 'solve_ivp'
    y : float
        dependent variables. 
        y[0] -> mass
        y[1] -> pressure
    """
    return y[1]/central_pressure - 1e-6