import numpy as np
import matplotlib.pyplot as plt
from .utils import *
from .eos import PressureDensityPolytropic, PressureDensityPiecewise

# get the eos from the class PressureDensityPiecewise (which uses PressureDensityPolytropic)

def get_dense_eos(key):
    """
    Parameters:
    ===========
    key : str
        string that identifies the tabulated eos
    -----------
    dense_eos : PressureDensityPiecewise
        eos of the core
    """
    # read 'dense' eos table for parameters
    # parameters[0] = log_10(pressure)
    # parameters[i] = gamma_i; for i = 1,2,3
    parameters = eos_lib[ key ]

    # dense eos starting pressure (in geom)
    p0 = (10**parameters[0])*conversion_dict['cgs']['pressure']['geom']

    #polytrope indices
    gamma1 = parameters[1] 
    gamma2 = parameters[2]
    gamma3 = parameters[3]

    #transition densities -> tabulated in the article
    trans1 = 2.7e14*conversion_dict['cgs']['density']['geom']
    trans2 = (10**14.7)*conversion_dict['cgs']['density']['geom']
    trans3 = (10**15.0)*conversion_dict['cgs']['density']['geom']

    #polytric constants -> formula from article of Read (eqs 6)
    k1 = p0/(trans2**gamma1)
    k2 = k1 * trans2**(gamma1-gamma2)
    k3 = k2 * trans3**(gamma2-gamma3)

    # we initialize a Polytropic Eos for each ranges.
    tropes = [PressureDensityPolytropic(k1, gamma1),
              PressureDensityPolytropic(k2, gamma2),
              PressureDensityPolytropic(k3, gamma3) ]
    
    # we pack transition densities
    transitions = [trans1, trans2, trans3]

    # we pack them in a Polytropic Piecewise class instance
    dense_eos = PressureDensityPiecewise(tropes, transitions)

    return dense_eos

# We now define the eos for the core, which is a SLy eos
# It's a particular instance of the class PressureDensityPiecewise
def sly_eos(): 

    # in cgs
    k_SLy = [6.80110e-9, 1.06186e-6, 5.32697e1, 3.99874e-8] # polytropic constants
    gamma_SLy = [1.58425, 1.28733, 0.62223, 1.35692] # polytropic indices
    transitions_SLy = [1.e3, 2.44034e7, 3.78358e11, 2.62780e12 ] # transition densities

    # we convert in geom
    k_SLy = [ (k_SLy[i] * conversion_dict['cgs']['lenght']['m']**(3*gamma_SLy[i] -1) 
                        * conversion_dict['cgs']['time']['geom']**(-2)
                        * conversion_dict['cgs']['mass']['geom'] ** (1-gamma_SLy [i]) )
                        for i in range(len(k_SLy)) ]
    transitions_SLy = [ x*conversion_dict['cgs']['density']['geom'] for x in transitions_SLy ]

    tropes = []
    trans = []

    prev_trope = None

    for (k, gamma, transition) in zip(k_SLy, gamma_SLy, transitions_SLy):
        #trope = PressureDensityPolytropic(k, gamma)
        trope =  PressureDensityPolytropic(k*C_CGS**2, gamma)
        tropes.append(trope)

        #correct transition depths to avoid jumps
        if not(prev_trope == None):
            rho_trans = (trope.k / prev_trope.k )**( 1.0/( prev_trope.gamma - trope.gamma ) )
        else:
            rho_trans = transition
        
        prev_trope = trope 
        trans.append(rho_trans)

    #Create crust using polytrope class
    SLyCrust = PressureDensityPiecewise(tropes, trans)
    return SLyCrust

# Smoothly glue core (dense eos) to SLy crust

def glue_crust_core_eos(key):
    """
    PARAMETERS:
    ===========
    key : str
        string that defines the dense_eos
    -----------
    eos : PressureDensityPiecewise 
        returned eos
    NOTES
    core : PressureDensityPiecewise 
        dense eos returnes by get_dense_eos
    crust : PressureDensityPiecewise 
        sly crust returned by sly_eos
    """
    core = get_dense_eos(key)
    crust = sly_eos()
    
    #unpack crust and core
    tropes_crust = crust.tropes 
    trans_crust  = crust.transitions

    tropes_core = core.tropes
    trans_core  = core.transitions

    #find transition depth
    rho_tr = (tropes_core[0].k / tropes_crust[-1].k )**( 1.0/( tropes_crust[-1].gamma - tropes_core[0].gamma ) )
        
    trans_core[0] = rho_tr

    #repack
    tropes = tropes_crust + tropes_core
    trans  = trans_crust  + trans_core

    for trope in tropes:
        trope.a = 0.0

    eos = PressureDensityPiecewise(tropes,trans)
    return eos

