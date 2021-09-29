import numpy as np
import matplotlib.pyplot as plt
import os
import compactobjects as cobj
from compactobjects import conversion_dict, eos_lib
#import palettable as pal
import math

# to style plots
#cmap = pal.cmocean.sequential.Matter_8.mpl_colormap 
colors = plt.rcParams["axes.prop_cycle"]()

# useful constants 

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

# create the folder where save the results

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, 'ReadOutput/')
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)

######################################
#       MAIN FOR NEUTRON STARS       #
######################################

# Study of Read EOS'

#  Build only a star with a tabulated value of central density

eos = cobj.glue_crust_core_eos('SLy')

ns = cobj.CompactStar(eos, value_type = 'Density')

rho0 = 3.5e15 * conversion_dict['cgs']['density']['geom'] 

r_tov, m_tov, p_tov = ns.structure_solver('TOV', rho0)

R_tov = r_tov[-1]
M_tov = m_tov[-1]

print('=========================================================')
print('SLy Neutron Star.')
print('Central density : ', rho0*conversion_dict['geom']['density']['cgs'], " g/cm^3")
# rho0, " g/cm^3")
print('---------------------------------------------------------')
print(' M_tov = ', M_tov, ' R_tov = ', R_tov)
print('=========================================================')

# create a figure of m(r) and P(r)

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.text(0.5, 1.07, "P(r) & m(r) of a SLy Neutron Star",
         horizontalalignment='center',
         fontsize=12,
         transform = ax.transAxes)
ax.plot(r_tov, p_tov, color="black", linestyle="-", linewidth=2,  label = 'P TOV')
ax.set_xlabel('r [km]',fontsize=10)
ax.set_ylabel(r'P [$dyne/cm^2$]', fontsize=10)
ax.minorticks_on()

ax2 = ax.twinx()
ax2.plot(r_tov, m_tov, color="black", linestyle="-.", label = 'm TOV')
ax2.plot(R_tov, M_tov, marker = 'o', color='grey', label='NS TOV mass')
ax2.set_ylabel(r"m [$M_{\odot}$]",fontsize=10)
ax2.minorticks_on()

fig.legend(loc="upper center", bbox_to_anchor=(0.5,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'sly_mp-vs-r.png',
            format='png',
            dpi=1000)

#  Build a sequence of 200 stars.

densities = np.logspace(13.5, 16.5, 300)*conversion_dict['cgs']['density']['geom']

# plot eos

eden = []
pressure = []

eden = [eos.eden_from_density(rho)*conversion_dict['geom']['pressure']['cgs'] for rho in densities]
pressure = [eos.pressure_from_density(rho)*conversion_dict['geom']['pressure']['cgs'] for rho in densities]

fig,ax = plt.subplots()
plt.title("SLy EOS")
ax.plot(pressure, eden, 
    color="black", linestyle=":", linewidth=2, label = r'$\epsilon$')
ax.set_xscale('log')
ax.minorticks_on()
ax.set_xlabel(r'P [$dyne/cm^2$]',fontsize=14)
ax.set_ylabel(r"$\epsilon$ [$erg/cm^3$]", fontsize=14)

fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'eos.png',
            format='png',
            dpi=1000) 

# solving structure eqs

R_star_tov, M_star_tov = ns.mass_vs_radius('TOV', densities)

fig,ax = plt.subplots()
plt.title("Mass-Radius of a SLy Neutron Star")
plt.rc('font', family='monospace')
ax.plot(R_star_tov, M_star_tov, color="black", linestyle=":", linewidth=2)
ax.set_xlabel('R [km]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)
ax.minorticks_on()
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'sly_mass-vs-radius.png',
            format='png',
            dpi=1000)

# plot Mass/Radius - central density 

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass/Radius vs Central Pressure in a SLy NS")
ax.set_xscale('log')
ax.minorticks_on()
ax.plot(densities*conversion_dict['geom']['density']['cgs'], M_star_tov, color="black", linestyle="-.", linewidth=2,  label = 'M-TOV')
ax.set_xlabel('p0 [$dyne/cm^2$]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)

ax2 = ax.twinx()
ax2.set_xscale('log')
ax2.minorticks_on()
ax2.plot(densities*conversion_dict['geom']['density']['cgs'], R_star_tov, color="black", linestyle=":", linewidth=2,  label = 'R-TOV')
ax2.set_ylabel(r"R [$km$]",fontsize=14)

fig.legend(loc='center right', bbox_to_anchor=(1,0.5), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'sly_mr-vs-p0.png',
            format='png',
            dpi=1000)

########################

# We now produce a comparison M-R plot between all the eos contained in eos_lib (contained in src/compactobjects/read_eos.py)

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Comparison between different EoS")
#ax.set_xlim(6, 25.0)
#ax.set_ylim(0.0, 3.0)
ax.set_xlabel(r'Radius $R$ [km]')
ax.set_ylabel(r'Mass $M$ [$M_{\odot}$]')
ax.minorticks_on()

i = 0
for key, value in eos_lib.items():
    print("Solving TOV equations for ", key, " EoS")
    eos = cobj.glue_crust_core_eos(key)
    ns = cobj.CompactStar(eos, value_type = 'Density')
    R_star_tov, M_star_tov = ns.mass_vs_radius('TOV', densities)
    
    # plot
    linestyle='solid'
    c = next(colors)["color"]

    ax.plot(R_star_tov, M_star_tov, color=c, linestyle=linestyle, alpha = 0.9, label = key)
    i += 1

fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
plt.savefig(results_dir+'mass-radius-comparison.pdf')