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