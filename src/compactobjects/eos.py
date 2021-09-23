import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from .constants import *
from .bisection import bisection

class ImplicitEos (object):

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

        self.x_range = np.linspace(-1,1,int(10000)) # default
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
        e_column = self.energy_fermi_momentum(x_range)
        p_column = self.pressure_fermi_momentum(x_range)
        return p_column, e_column

    def set_x_range (self, new_x_range):
        """
        Redefine self_x_range
        Parameters:
        ===========
        new_x_range : array of float
            for real high pressure we need real high 
            fermi momenta (x ~ 1e2/1e3). The higher 
            number of points, the higher precision.
            The range must be symmetric w.r.t. 0!
        """
        self.x_range = new_x_range
        
    """
    # version used by the article: it works as cubicspline but much slower
    # to use this, replace lines 47->76 with the following function and uncomment line 5

    def eden_from_pressure (self, pressure): 
        pressure_shifted = lambda x : self.pressure_fermi_momentum(x) - pressure
        x_p = bisection(pressure_shifted, -100, 100, 10000)  # really slow
        eden = self.energy_fermi_momentum(x_p)
        return eden
    """
    def plot(self, case = 'Rel'):
        """
        Function that plots the points p and e computed from the explicit function e(x) and p(x).
        We generate x between 0 and 0.1 for the non-relativistic case, 0-1 for the relativistic one and between 1 
        and 10 for the ultrarelativistic case.

        Parameters
        ============
        case : string
            choose between the three cases
        """
        fig, ax = plt.subplots()
        plt.title("Energy density - Pressure EOS")
        if case == 'NR' : x_range = np.linspace(0,0.1,10000)
        elif case == 'Rel': x_range = np.linspace(0,1,10000)
        elif case == 'UR': x_range = np.linspace(1,10,10000)
        else : 
            raise ValueError("'case' can be only 'NR', 'Rel' or 'UR' depending if you want to plot the EOS" 
                " in the non-relativistic, full relativistic or ultra-relativistic case respectively.")
        p_range = self.pressure_fermi_momentum(x_range)
        e_range = self.energy_fermi_momentum(x_range)
        ax.plot(p_range*conversion_dict['geom']['pressure']['cgs'], e_range*conversion_dict['geom']['pressure']['cgs'], 
            color="black", linestyle="-." ) 
        ax.set_xlabel(r'P [$dyne/cm^2$]', fontsize=14)
        ax.set_ylabel(r'$\epsilon$  [$dyne/cm^2$]',fontsize=14)

        fig.savefig('eos.png',
            format='png',
            dpi=1000,
            bbox_inches='tight')    

class PolytropicEos(object): 

    def __init__ (self, k, gamma):
        """
        Given a value k and a polytropic index gamma return energy density in function of pressure.
        Needs geometric units as input (sse 'conversion_dict' taken from constants.py. See there the usage of the dictionary.)
        Returns in CGS units.

        Parameters
        ==========
        k: float
            tropic coefficiente
        gamma: float
            adiabatic index
        n : float
            polytropic index
        """
        self.type = 'PolytropicEos'
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

    def plot(self, plot_range):
        """
        Function that plots the points p and e computed from the explicit function e(x) and p(x).
        We generate x between 0 and 0.1 for the non-relativistic case, 0-1 for the relativistic 
        one and in [1,10] for the ultrarelativistic case.
        Use 'conversion_dict' taken from constants.py. See there the usage of the dictionary.

        Parameters
        ============
        plot_range : array of float
            array of pressures that we want plot
        """
        fig, ax = plt.subplots()
        plt.title("Energy density - Pressure EOS")
        p_range = plot_range*conversion_dict['geom']['pressure']['geom']
        e_range = self.eden_from_pressure(p_range)
        ax.plot(p_range*conversion_dict['geom']['pressure']['cgs'], e_range*conversion_dict['geom']['pressure']['cgs'], 
            color="black", linestyle="-." ) 
        ax.set_xlabel(r'P [$dyne/cm^2$]', fontsize=14)
        ax.set_ylabel(r'$\epsilon$  [$dyne/cm^2$]',fontsize=14)

        fig.savefig('eos.png',
            format='png',
            dpi=1000,
            bbox_inches='tight')    