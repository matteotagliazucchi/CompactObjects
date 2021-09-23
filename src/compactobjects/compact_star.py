import numpy as np 
from .runge_kutta import rk_odesolver
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
        self.radii = np.linspace(1e-6,1e5,10000) # in m (geom)
        self.eos = eos

    def change_radii_range (self, max_radius):
        """
        Function that allows to change the max radius of self.radii

        Parameters:
        ===========
        max_radius : array of floats
            new max radius. min one is by default 1e-6 m 
        """
        self.radii = np.linspace(1e-6, max_radius, 100000)

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
        Utility to set initial condition.

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

    def structure_solver(self, eqs_type, central_pressure):
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
        #print("Solving ", eqs_type, " equations...")
        
        if self.eos.type == 'ImplicitEos' : 
            # we change the default fermi energy momenta range for eos in order to surely include all the pressures
            pressure_shifted = lambda x : self.eos.pressure_fermi_momentum(x) - central_pressure
            x_p = bisection(pressure_shifted, -1e3, 1e3, int(1e5))
            x_new = np.linspace(-int(x_p)-1, int(x_p)+1, 10000)
            self.eos.set_x_range(x_new)
        
        initial = self.set_initial_conditions(central_pressure)

        if (eqs_type == 'Newton'):

            solutions, check_event = rk_odesolver(self.newton_eqs, initial, self.radii, stop_event = self.found_radius)
            
            if check_event == False :
                # print ("Findinding where the pressure reaches the default minium value...")               
                # if the pressure has not reached a negative value we search when it goes under a default value w.r.t the central pressure
                for i, p in enumerate(solutions[:,1]):
                    if ( (p/central_pressure) < 1.e-6 ) :
                        m_out = solutions[:i,0]*conversion_dict['geom']['mass']['m_sol']
                        p_out = solutions[:i,1]*conversion_dict['geom']['pressure']['cgs']
                        r_out = self.radii[:len(p_out)]*conversion_dict['geom']['lenght']['km'] 
                        return r_out, m_out, p_out
                
                raise ValueError("Minimum pressure not reached. Unable to determine star mass and radius.")
                quit()
            
            else : # we reached a minimum pressure 

                # we take the second to last element because the last element is an immaginary pressure (and a NaN mass)!
                m_out = solutions[:-2,0]*conversion_dict['geom']['mass']['m_sol']  
                p_out = solutions[:-2,1]*conversion_dict['geom']['pressure']['cgs']
                r_out = self.radii[:len(p_out)]*conversion_dict['geom']['lenght']['km'] 
                return r_out, m_out, p_out
        
        elif (eqs_type == 'TOV'):
            
            solutions, check_event = rk_odesolver(self.tov_eqs, initial, self.radii, stop_event = self.found_radius)
            
            if check_event == False :
                #print ("Findinding where the pressure reaches the default minium value...")      
                # if the pressure has not reached a negative value we search when it goes under a default value
                for i, p in enumerate(solutions[:,1]):
                    if ( (p/central_pressure) < 1.e-6 ) :
                        m_out = solutions[:i,0]*conversion_dict['geom']['mass']['m_sol']
                        p_out = solutions[:i,1]*conversion_dict['geom']['pressure']['cgs']
                        r_out = self.radii[:len(p_out)]*conversion_dict['geom']['lenght']['km'] 
                        return r_out, m_out, p_out
                
                raise ValueError("Minimum pressure not reached. Unable to determine star mass and radius.")
                quit()
            
            else : # we reached a minimum pressure 

                # we take the second to last element because the last element is an immaginary pressure (and a NaN mass)!
                m_out = solutions[:-2,0]*conversion_dict['geom']['mass']['m_sol'] 
                p_out = solutions[:-2,1]*conversion_dict['geom']['pressure']['cgs']
                r_out = self.radii[:len(p_out)]*conversion_dict['geom']['lenght']['km'] 
                return r_out, m_out, p_out
        else:
            raise ValueError("Must choose 'Newton' or 'TOV' in 'eqs_type'")
            quit()

    def found_radius(self,t,y):
        """
        Search when the pressure becomes a NaN. 
        Useful for polytropic eos, because in this case the pressure out of the star surface becomes immaginary (numerically a "NaN")

        Parameters:
        ===========
        t : float
            independent variable. 
            Not used, but necessary to have the right signature required by 'rk_odesolver'
        y : float
            dependent variables. 
            y[0] -> mass
            y[1] -> pressure
        """
        return (np.isnan(y[1]))

    def mass_vs_radius (self, eqs_type, central_pressures):
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
            r_star , m_star, p_star = self.structure_solver(eqs_type, pc)
            R_star[i] = r_star[-1]
            M_star[i] = m_star[-1]

        return R_star, M_star


