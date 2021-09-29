import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from .utils import *
from .bisection import bisection

###################
#   IMPLICIT EOS #
##################
class ImplicitEos (object):

    x_range = np.linspace(-1,1,int(10000))

    def __init__(self, energy_fermi_momentum, pressure_fermi_momentum):

        """
        Class that, given the explicit expression of the energy and the pressure 
        in function of x = k_f/(m*c), where k_f is the Fermi momentum,
        returns the interpolated eos, i.e the interpolated 'e = e(p)'
        Needs geometric units as input (use 'conversion_dict' taken from constants.py. 
        See there the usage of the dictionary). Returns in CGS units.

        Parameters
        ===========
        energy_fermi_momentum : callable -> can be also a lambda function
            energy density in function of x = k_f/(m*c) 
        pressure_fermi_momentum : callable -> can be also a lambda function
            pressure in function of x = k_f/(m*c) 
        x_range : array of float
            points x = k_f/(m*c) where to compute e(x) and p(x)
            default value is [-1,1] (and 10000 points) works for quite 
            large pressure, but for real high pressure we need real high 
            fermi momenta (x ~ 1e2/1e3). The higher number of points, the 
            higher precision.
        p_column : array of float
            pressure computed from the full expression for x in x_range
        e_column : array of float
            energy computed from the full expression for x in x_range
        eden_from_pressure : callable
            interpolation between p_column and e_column
            Needs a pressure input where to compute 'eden'
            Basic usage: 'eden = eos.eden_from_pressure(pressure)'
        -----------
        NOTES:  we generates the Fermi Momnetum between -1 and 1 to make 
                CubicSpline works correctly, indeedin 'rk_fixedstep' 
                we need negative pressures sometimes to computed the step 
                and so CubicSline needs x < 0.
        """
        self.type = 'ImplicitEos'
        self.energy_fermi_momentum = energy_fermi_momentum
        self.pressure_fermi_momentum = pressure_fermi_momentum

        self.p_column, self.e_column = self.eden_pressure_table (self.x_range)

        self.eden_from_pressure = CubicSpline(self.p_column, self.e_column)

    def eden_pressure_table (self, x_range):
        """
        Build p_column and e_column given a range of Fermi momenta 

        Parameters:
        ===========
        x_range : array of float
            array of Fermi momenta where compute p(x) and e(x) 
        """
        print('Costruisco tabelle')
        e_column = self.energy_fermi_momentum(x_range)
        p_column = self.pressure_fermi_momentum(x_range)
        return p_column, e_column
    
    """
    # version used by the article: it works as cubicspline but much slower
    # to use this, replace lines 47->76 with the following function and uncomment line 5

    def eden_from_pressure (self, pressure): 
        pressure_shifted = lambda x : self.pressure_fermi_momentum(x) - pressure
        x_p = bisection(pressure_shifted, -100, 100, 10000)  # really slow
        eden = self.energy_fermi_momentum(x_p)
        return eden
    """  

#################################
#  PRESSUE-EDEN POLYTROPIC EOS  #
#################################

class PressureEdenPolytropic(object): 

    def __init__ (self, k, gamma):
        """
        Given the Polytropic eos P = K*\epsilon^\gamma, that relates the pressure P and the the enegy density \epsilon.
        Return energy density in function of pressure
        Needs geometric units as input (sse 'conversion_dict' taken from constants.py. See there the usage of the dictionary.)
        Returns in CGS units.

        Parameters
        ==========
        type : str
            to distinguish among different classes
        k: float
            tropic coefficiente
        gamma: float
            adiabatic index
        """
        self.type = 'PressureEdenPolytropic'
        self.k = k
        self.gamma = gamma   

    def eden_from_pressure(self, pressure):
        """
        Polytropic  eden = (P/k)**(1/gamma)
        
        Parameters:
        ===========
        pressure: float
            pressure 
        -----------
        eden : float
            energy density
        """
        eden = (pressure/self.k)**(1/self.gamma)
        return eden

####################################
#  PRESSUE-DENSITY POLYTROPIC EOS  #
####################################

class PressureDensityPolytropic (object):
    
    #transition continuity constant 
    a = 0.0
    
    def __init__ (self, k, gamma):
        """
        Polytropic eos P = K*\rho^\gamma, that relates the pressure P and the the rest mass density \rho.
        Needs geometric units as input (sse 'conversion_dict' taken from constants.py. See there the usage of the dictionary.)
        Returns in CGS units.

        Parameters
        ==========
        type : str
            to distinguish among different classes
        k: float
            tropic coefficiente
        gamma: float
            adiabatic index
        n : float
            polytropic index
        """        
        self.type = 'PressureDensityPolytropic'
        #self.k = k/C_CGS**2
        self.k = k
        self.gamma = gamma  
        self.n = 1/(self.gamma - 1) 

    def pressure_from_density (self, density):
        """
        Polytropic P = k*rho**(gamma)
        
        Parameters:
        ===========
        density : float
            rest mass density
        -----------
        pressure: float
            pressure 
        """
        pressure = self.k*(density**self.gamma)
        #pressure = (C_CGS**2)*self.k*(density**self.gamma)
        return pressure  

    def density_from_pressure(self, pressure):
        """
        Polytropic rho = (P/k)**(1/gamma)
        
        Parameters:
        ===========
        pressure: float
            pressure 
        -----------
        rho : float
            rest mass density
        """
        if pressure < 0.0:
            return 0.0
        #rho = (pressure/C_CGS**2/self.k)**(1/self.gamma)
        rho = (pressure/self.k)**(1/self.gamma)
        return rho

    def eden_from_density (self, density):
        """
        eden = (1+a)*rho + n*k*rho^\gamma
        
        Parameters:
        ===========
        density: float
            rest mass density 
        -----------
        eden : float
            energy density
        """ 
        eden = (1.+self.a)*density + self.n*self.k*(density**self.gamma)
        # maybe a problem here due to a factor c**2? check article
        return eden

    def eden_from_pressure(self, pressure):
        """
        Polytropic  eden = (P/k)**(1/gamma)
        
        Parameters:
        ===========
        pressure: float
            pressure 
        -----------
        eden : float
            energy density
        """
        density = self.density_from_pressure(pressure) 
        eden = self.eden_from_density(density)
        return eden
#-----------------------------------------------------------------------------------------------------------#
class PressureDensityPiecewise (object):

    def __init__ (self, tropes, trans, prev_trope = None ):
        """
        Piecewise Polytropic eos P = K_i*\rho^\gamma_i, that relates the pressure P and the the rest mass density \rho.
        Needs geometric units as input (sse 'conversion_dict' taken from constants.py. See there the usage of the dictionary.)
        Returns in CGS units.

        Parameters
        ==========
        type : str
            to distinguish among different classes
        tropes : PressureDensityPolytropic
            different pressure-rho polytropic eos
        trans: aray of float
            densities where the we have to change the polytropic
        prev_trope : PressureDensityPolytropic
            pressure-density polytropic of the previous range
        pressures : array of float
            pressures corresponding to the densities transition
        edens : array of float
            eden corresponding to the densities transition
        """
        self.type = 'PressureDensityPiecewise'
        self.tropes      = tropes
        self.transitions = trans

        self.pressures  = []
        self.edens  = []

        for (trope, transition) in zip(self.tropes, self.transitions):

                if not( prev_trope == None ):
                    trope.a = self.set_trans_const( prev_trope, trope, transition )
                else:
                    transition = 0.0

                eden = trope.eden_from_density(transition) 
                pressure = trope.pressure_from_density(transition)

                # eden and pressures corresponding to the densities transition values

                self.pressures.append( pressure )
                self.edens.append( eden )

                prev_ed = eden
                prev_tr = transition
                prev_trope = trope

    def set_trans_const(self, prev_trope, trope, transition):
        """
        Parameters:
        ===========
        trope : PressureDensityPolytropic
            different pressure-rho polytropic eos
        transition: aray of float
            densities where the we have to change the polytropic
        prev_trope : PressureDensityPolytropic
            pressure-density polytropic of the previous range
        """
        return prev_trope.a + (prev_trope.k/(prev_trope.gamma - 1))*transition**(prev_trope.gamma-1) - (trope.k/(trope.gamma - 1))*transition**(trope.gamma-1)

    def find_interval_given_density(self, density):
        """
        Parameters:
        ===========
        density : float
            density
        -----------
        trope : PressureDensityPolytropic
            polytropic press-density eos depending on the density range where we are
        """

        if density <= self.transitions[0]:
            return self.tropes[0]

        for x in range( len(self.transitions) - 1 ):
            if self.transitions[x] <= density < self.transitions[x+1]:
                return self.tropes[x]

        return self.tropes[-1]

    def find_interval_given_pressure(self, pressure):
        """
        Parameters:
        ===========
        pressure : float
            pressure
        -----------
        trope : PressureDensityPolytropic
            polytropic press-density eos depending on the pressure range where we are
        """
        if pressure <= self.pressures[0]:
            return self.tropes[0]

        for x in range( len(self.pressures) - 1):
            if self.pressures[x] <= pressure < self.pressures[x+1]:
                return self.tropes[x]

        return self.tropes[-1]

    def density_from_pressure(self, pressure):
        """
        Parameters:
        ===========
        pressure: float
            pressure 
        -----------
        density : float
            rest mass density
        """
        trope = self.find_interval_given_pressure(pressure)
        return trope.density_from_pressure(pressure)

    def pressure_from_density(self, density):
        """
        Parameters:
        ===========
        density : float
            rest mass density
        -----------
        pressure: float
            pressure 
        """
        trope = self.find_interval_given_density(density)
        return trope.pressure_from_density(density)

    def eden_from_density(self, density):
        """
        Parameters:
        ===========
        density : float
            rest mass density
        -----------
        eden: float
            energy density  
        """
        trope = self.find_interval_given_density(density)
        return trope.eden_from_density(density)

    def eden_from_pressure(self, pressure):
        """
        Parameters:
        ===========
        pressure : float
            pressure
        -----------
        eden: float
            energy density  
        """
        trope = self.find_interval_given_pressure(pressure)
        return trope.eden_from_pressure(pressure)

    #vectorized version
    def pressures_from_densities(self, densities):
        press = []
        for rho in densities:
            press.append( self.pressure_from_density(rho) )
        return press