from scipy.constants import speed_of_light, gravitational_constant, pi, hbar, h, m_n, m_e
import astropy.constants as const

solar_mass = const.M_sun.value  # kg

C_SI = speed_of_light  # m s^-1
C_CGS = C_SI*100 # cm s^-1
G_SI = gravitational_constant  # m^3 kg^-1 s^-2
G_CGS = G_SI*1000 # cm^3 g^-1 s^-2
MSUN_SI = solar_mass  # kg
MSUN_CGS = solar_mass*1000  # g
HBAR_SI = hbar # J s
HBAR_CGS = HBAR_SI*1e7 # ergs s

# Stores conversions from geometerized to CGS or SI units systems and vice versa.
# Conversions taken from Wald book.
# Geometrized units are defined as c = G = 1 and we arbitrarly choose to put c and G in SI units so that we have 1 m = 1 s and so on.
# In geom units everything is expressed in power of m. Again, see Wald General relativity book.
# usage: quantity_converted = initial_quantity*conversion_dict['initial units systme']['physical quantity']['final units system']

conversion_dict = { 'geom': {'pressure': {'cgs': 10*(C_SI ** 4.) / G_SI,
                                          'si': (C_SI ** 4.) / G_SI,
                                          'geom': 1.
                                          },

                             'energy_density': {'cgs': 10*(C_SI ** 4.) / G_SI,
                                                'si': (C_SI ** 4.) / G_SI,
                                                'geom': 1.
                                                },

                             'density': {'cgs': 1e-3*(C_SI ** 2.) / G_SI,
                                         'si': (C_SI ** 2.) / G_SI,
                                         'geom': 1.
                                         },

                             'mass': {'g': 1e3*(C_SI ** 2.) / G_SI,
                                      'kg': (C_SI ** 2.) / G_SI,
                                      'geom': 1.,
                                      'm_sol': (C_SI ** 2.) / G_SI / MSUN_SI
                                      },

                             'lenght': {'cm': 100.,
                                        'm': 1,
                                        'km': 1e-3
                                        },

                             'time': {'cgs': 1/C_SI,
                                      'si': 1/C_SI,
                                      'geom': 1.
                                      },

                             'energy': {'cgs': 1e7*(C_SI ** 4.) / G_SI, 
                                        'si': (C_SI ** 4.) / G_SI, 
                                        'geom': 1.
                                        }
                            }, # end of conversion from GEOM

                    'cgs' : { 'pressure': {'cgs': 1,
                                         'si': 0.1 ,
                                         'geom': 0.1*G_SI/(C_SI ** 4.)
                                        },

                             'energy_density': { 'cgs': 1,
                                                'si': 0.1,
                                                'geom': 0.1*G_SI/(C_SI ** 4.)
                                                },

                             'density': { 'cgs': 1.,
                                         'si': 1e3,
                                         'geom': 1e3*G_SI/(C_SI ** 2.)
                                         },

                             'mass': { 'g': 1. ,
                                      'kg': 1e-3,
                                      'geom': 1e-3*G_SI/(C_SI ** 2.) ,
                                      'm_sol': 1/MSUN_CGS
                                      },

                             'lenght': { 'cm': 1.,
                                        'm': 1e-2,
                                        'km': 1e-5
                                        },

                             'time': { 'cgs': 1,
                                      'si': 1,
                                      'geom': C_SI
                                      },

                             'energy': { 'cgs': 1., 
                                        'si': 1e-7, 
                                        'geom': 1e-7*G_SI/(C_SI ** 4.)
                                        }
                            }, # end of conversion from CGS

                    'si' : {'pressure': {'cgs': 10.,
                                          'si': 1.,
                                          'geom': G_SI/(C_SI ** 4.)
                                          },

                             'energy_density': {'cgs': 10.,
                                                'si':1.,
                                                'geom': G_SI/(C_SI ** 4.)
                                                },

                             'density': {'cgs': 1e-3,
                                         'si':1.,
                                         'geom': (C_SI ** 2.) / G_SI,
                                         },

                             'mass': {'g': 1000.,
                                      'kg': 1.,
                                      'geom':  G_SI/(C_SI ** 2.),
                                      'm_sol': 1 / MSUN_SI
                                      },

                             'lenght': {'cm': 100.,
                                        'm': 1.,
                                        'km': 1e-3
                                        },

                             'time': {'cgs': 1.,
                                      'si': 1.,
                                      'geom': C_SI,
                                      },

                             'energy': {'cgs': 1e7, 
                                        'si': 1., 
                                        'geom':  G_SI/(C_SI ** 4.)
                                        }
                             } # end of conversion from SI

                    } # end of dictionary


########################
#  READ TABULATED EOS  #
########################

#Dictionary from Read et al 2009 
# all M_max < 2Msun commented out
# units in cgs!
eos_lib = {
#    'PAL6'  :[ 34.380,  2.227,  2.189,  2.159, 'npem' ],
    'SLy'   :[ 34.384,  3.005,  2.988,  2.851, 'npem' ],
#    'APR1'  :[ 33.943,  2.442,  3.256,  2.908, 'npem' ],
#    'APR2'  :[ 34.126,  2.643,  3.014,  2.945, 'npem' ],
    'APR3'  :[ 34.392,  3.166,  3.573,  3.281, 'npem' ],
    'APR4'  :[ 34.269,  2.830,  3.445,  3.348, 'npem' ],
#    'FPS'   :[ 34.283,  2.985,  2.863,  2.600, 'npem' ],
    'WFF1'  :[ 34.031,  2.519,  3.791,  3.660, 'npem' ],
    'WFF2'  :[ 34.233,  2.888,  3.475,  3.517, 'npem' ],
#    'WFF3'  :[ 34.283,  3.329,  2.952,  2.589, 'npem' ],
#    'BBB2'  :[ 34.331,  3.418,  2.835,  2.832, 'npem' ],
#    'BPAL12':[ 34.358,  2.209,  2.201,  2.176, 'npem' ],
    'ENG'   :[ 34.437,  3.514,  3.130,  3.168, 'npem' ],
    'MPA1'  :[ 34.495,  3.446,  3.572,  2.887, 'npem' ],
    'MS1'   :[ 34.858,  3.224,  3.033,  1.325, 'npem' ],
#    'MS2'   :[ 34.605,  2.447,  2.184,  1.855, 'npem' ],
    'MS1b'  :[ 34.855,  3.456,  3.011,  1.425, 'npem' ],
#    'PS'    :[ 34.671,  2.216,  1.640,  2.365, 'meson' ],
#    'GS1a'  :[ 34.504,  2.350,  1.267,  2.421, 'meson' ],
#    'GS2a'  :[ 34.642,  2.519,  1.571,  2.314, 'meson' ],
#    'BGN1H1':[ 34.623,  3.258,  1.472,  2.464, 'hyperon' ],
#    'GNH3'  :[ 34.648,  2.664,  2.194,  2.304, 'hyperon' ],
#    'H1'    :[ 34.564,  2.595,  1.845,  1.897, 'hyperon' ],
#    'H2'    :[ 34.617,  2.775,  1.855,  1.858, 'hyperon' ],
#    'H3'    :[ 34.646,  2.787,  1.951,  1.901, 'hyperon' ],
    'H4'    :[ 34.669,  2.909,  2.246,  2.144, 'hyperon' ],
#    'H5'    :[ 34.609,  2.793,  1.974,  1.915, 'hyperon' ],
#    'H6a'   :[ 34.593,  2.637,  2.121,  2.064, 'hyperon' ],
#    'H7'    :[ 34.559,  2.621,  2.048,  2.006, 'hyperon' ],
#    'PCL2'  :[ 34.507,  2.554,  1.880,  1.977, 'hyperon' ],
#    'ALF1'  :[ 34.055,  2.013,  3.389,  2.033, 'quark' ],
    'ALF2'  :[ 34.616,  4.070,  2.411,  1.890, 'quark' ],
#    'ALF3'  :[ 34.283,  2.883,  2.653,  1.952, 'quark' ],
#    'ALF4'  :[ 34.314,  3.009,  3.438,  1.803, 'quark' ]
}
