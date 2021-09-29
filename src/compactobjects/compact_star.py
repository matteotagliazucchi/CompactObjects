import numpy as np 
from scipy.integrate import solve_ivp
from .utils import  *

class CompactStar (object):

    def __init__ (self, eos, value_type = 'Pressure'):
        """
        Given an eos return mass and radius of the star (given also an external central pressure/density).

        Parameters
        ==========
        eos: callable function
            polytropicEOS() instance
        value_type : str
            'Pressure' or 'Density' depending id we provide a central pressure or a central density
        radii = array of float
            range of radii where to search the real radius of the star via tov eqs.
            this range must contain the expected radius of the star: for a ns we must have 
            25 km (25e3 m) as maximum radius, while for a white dwarf we shold have 5e7 km.
            range expressed in geom units, i.e in meters. Default value set to work with Neutron Stars.
        """       
        self.radii = np.arange(1e-3,25e3,1) # in m (geom); default value
        self.eos = eos
        self.set_integration_options()
        if value_type == 'Pressure': 
            self.value_type = value_type 
        elif value_type == 'Density':
            if self.eos.type == 'ImplicitEos' or self.eos.type == 'PressureEdenPolytropic':
                raise ValueError('No support for density ad central value in ', self.eos.type, ' class.')
            self.value_type = value_type 
        else :
            raise ValueError('CompactStar works only if you provide a central value of Pressure or Density')

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
        #dy[1] = -G_CGS*(eden+p/C_CGS**2)*(m + 4*pi*r**3*p/C_CGS**2)/(r*(r-2*G_CGS*m/C_CGS**2))
        dy[1] = - (eden+p)*(m + 4*pi*r**3*p)/(r*(r-2*m))

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
  
    def set_initial_conditions(self, central_value): 
        """
        Utility to set initial conditions.

        Parameters:
        ===========
        central_value : float
            pressure or density at r = r_min
        -----------
        m0, p0: float
                initial mass and central pressure
        """
        if self.value_type == 'Density':
            rho0 = central_value
            p0 = self.eos.pressure_from_density(rho0)
        else: 
            p0 = central_value
        e0 = self.eos.eden_from_pressure(p0)
        m0 = 4.*pi*e0*self.radii[0]**3 
        return m0, p0

    def set_integration_options (self, method = 'RK45', rtol = 1.0e-4 , 
                atol = 0.0,  max_step = np.inf, dense_output = False):
        self.max_step = max_step
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.dense_output = dense_output

    def structure_solver(self, eqs_type, central_value):
        """
        Solve structure equations given a central pressure.

        Parameters:
        ===========
        central_value : float
            pressure or density at r = r_min
        eqs_type : string
            specify if we solve TOV or Newton eqs. Must be 'TOV' or 'Newton'
        max_step : float
            max step allowed for solve_ivp
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
            # too heavy in case of implicit eos to look for zero pressure, which never occurs
            stop = lambda r,y : self.found_radius_implicit(r, y, central_value) 
            stop.terminal = True  
        else: 
            stop = lambda r,y : self.found_radius(r, y)
            stop.terminal = True  

        initial = self.set_initial_conditions(central_value)
        
        r_span = [self.radii[0], self.radii[-1]]        
        
        solutions = solve_ivp(eqs_dict[eqs_type], r_span, initial, method = self.method , rtol = self.rtol,
                    atol=self.atol, max_step = self.max_step, dense_output = self.dense_output, events = stop)
        
        r_out = solutions.t*conversion_dict['geom']['lenght']['km'] 
        m_out = solutions.y[0]*conversion_dict['geom']['mass']['m_sol']  
        p_out = solutions.y[1] *conversion_dict['geom']['pressure']['cgs'] 
        
        if (solutions.status != 1): # pressure is not negative, so we check if it goes under a certain value
            for i,p in enumerate(p_out):
                if  p/initial[1] < 1e-6 :
                    r_out = r_out[:i]
                    m_out = m_out[:i]
                    p_out = p_out[:i]
                    return r_out, m_out, p_out

        return r_out, m_out, p_out

    def mass_vs_radius (self, eqs_type, central_values):
        """
        Returns arrays containing the Mass and the Radius of the compact stars, given a range of central pressures.
        A couple of values (R,M) for eache pressure

        Parameters:
        ===========
        central_values : array of float
            central pressures/densities
        eqs_type : string
            specify if we solve TOV or Newton eqs. Must be 'TOV' or 'Newton'
        -----------
        R_out, M_out : array of floats
            (R,M) of the star for each given central value. 
        """        
        R_star = np.empty_like(central_values)
        M_star = np.empty_like(central_values)

        for i, pc in enumerate(central_values):
            r_star , m_star, p_star = self.structure_solver(eqs_type, pc)
            R_star[i] = r_star[-1]
            M_star[i] = m_star[-1]

        return R_star, M_star

    def found_radius(self, t, y): # look here
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
        central_pressure : float
            central pressure
        """
        return y[1] # check when pressure is negative
            
    def found_radius_implicit(self, t, y, central_value):
        return y[1]/central_value - 1e-6


    # END OF CompactStar class
#===================================================================================#


