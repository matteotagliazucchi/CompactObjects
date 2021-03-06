import numpy as np 
import matplotlib.pyplot as plt
import os
import compactobjects as cobj
from compactobjects import conversion_dict

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
results_dir = os.path.join(script_dir, 'NSOutput/')
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)

######################################
#       MAIN FOR NEUTRON STARS       #
######################################

# POLYTROPIC NON RELATIVITISTIC CASE

#  Build only a star with a tabulated value of central pressure

gamma = 5/3
k = 6.483e-26*(conversion_dict['cgs']['lenght']['m']**2)*(conversion_dict['cgs']['energy']['geom'] ** (-2/3)) 
e0 = 1.603e38 * conversion_dict['cgs']['pressure']['geom'] 
p0 = e0 *1e-5

eos = cobj.PressureEdenPolytropic(k, gamma) 
ns = cobj.CompactStar(eos)

ns.set_radii_range(np.arange(1e-6, 1e5, 0.1))

r_newton , m_newton, p_newton = ns.structure_solver('Newton', p0)
R_newton = r_newton[-1]
M_newton= m_newton[-1]

r_tov , m_tov, p_tov = ns.structure_solver('TOV', p0)
R_tov = r_tov[-1]
M_tov = m_tov[-1]

print('=========================================================')
print('Pure Neutron Star in the non-relativistic polytropic case.')
print('Central pressure : ', p0*conversion_dict['geom']['pressure']['cgs'], " dyne/cm^2")
print('---------------------------------------------------------')
print(' M_newton = ', M_newton, ' R_newton = ', R_newton,
      '\n', ' M_tov = ', M_tov, ' R_tov = ', R_tov)
print('=========================================================')

# create a figure of m(r) and P(r)

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.text(0.5, 1.07, "P(r) & m(r) of a pure non relativistic neutron star",
         horizontalalignment='center',
         fontsize=12,
         transform = ax.transAxes)
ax.plot(r_newton, p_newton, color="blue", linestyle="-", linewidth=1, label = 'P Newton')
ax.plot(r_tov, p_tov, color="black", linestyle="-", linewidth=2,  label = 'P TOV')
ax.set_xlabel('r [km]',fontsize=14)
ax.set_ylabel(r'P [$dyne/cm^2$]', fontsize=14)
ax.minorticks_on()

ax2 = ax.twinx()
ax2.minorticks_on()
ax2.plot(r_newton, m_newton,color="blue", linestyle=":", label = 'm Newton')
ax2.plot(r_tov, m_tov, color="black", linestyle="-.", label = 'm TOV')
ax2.plot(R_newton, M_newton, marker = 'o', linestyle="", color='green', label='NS Newton mass')
ax2.plot(R_tov, M_tov, marker = 'o', color='red', linestyle="", label='NS TOV mass')
ax2.set_ylabel(r"m [$M_{\odot}$]",fontsize=14)

fig.legend(loc="upper center", bbox_to_anchor=(0.5,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'nrns_mp-vs-r.pdf',
            format='pdf',
            dpi=1000)

#  Build a sequence of 500 stars. 

p0_min = 1e30*conversion_dict['cgs']['pressure']['geom']
p0_max = 1e33*conversion_dict['cgs']['pressure']['geom']
pressures = np.linspace(p0_min, p0_max, 200)

R_star_tov, M_star_tov = ns.mass_vs_radius('TOV', pressures)
R_star_newton, M_star_newton = ns.mass_vs_radius('Newton', pressures)
# plot Mass-Radius

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass-Radius of a pure non relativistic NS")
ax.plot(R_star_newton, M_star_newton, color="blue", linestyle="-.", linewidth=1, label = 'Newton')
ax.plot(R_star_tov, M_star_tov, color="black", linestyle="-.", linewidth=2,  label = 'TOV')
ax.set_xlabel('R [km]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)
ax.minorticks_on()

fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'nrns_mass-vs-radius.pdf',
            format='pdf',
            dpi=1000)

# plot Mass/Radius - central pressure 
fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass/Radius vs Central Pressure in a pure non relativistic NS")
ax.set_xscale('log')
ax.minorticks_on()
ax.plot(pressures*conversion_dict['geom']['pressure']['cgs'], M_star_newton, color="blue", linestyle="-.", linewidth=1, label = 'M-Newton')
ax.plot(pressures*conversion_dict['geom']['pressure']['cgs'], M_star_tov, color="black", linestyle="-.", linewidth=2,  label = 'M-TOV')
ax.set_xlabel('p0 [$dyne/cm^2$]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)

ax2 = ax.twinx()
ax2.set_xscale('log')
ax2.minorticks_on()
ax2.plot(pressures*conversion_dict['geom']['pressure']['cgs'], R_star_newton, color="blue", linestyle=":", linewidth=1, label = 'R-Newton')
ax2.plot(pressures*conversion_dict['geom']['pressure']['cgs'], R_star_tov, color="black", linestyle=":", linewidth=2,  label = 'R-TOV')
ax2.set_ylabel(r"R [$km$]",fontsize=14)

fig.legend(loc='upper center', bbox_to_anchor=(0.5,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'nrns_mr-vs-p0.pdf',
            format='pdf',
            dpi=1000)