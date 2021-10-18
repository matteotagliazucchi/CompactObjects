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

# FULL RELATIVITISTIC CASE

#  Build only a star with a tabulated value of central pressure

e0 = (m_n**4)*(C_SI**5)
e0 = e0/((pi**2)*(HBAR_SI**3))
e0 = e0*conversion_dict['si']['energy_density']['geom']

e = lambda x :  e0*( (1/8)*((2*x**3 + x)*(1 + x**2)**(1/2) - np.arcsinh(x)) ) 
p = lambda x : (e0/24)*((2*x**3 - 3*x)*(1 + x**2)**(1/2) + 3*np.arcsinh(x))

eos = cobj.ImplicitEos(e,p)
ns = cobj.CompactStar(eos)

p0 = 3.5e35*conversion_dict['cgs']['pressure']['geom']

r_newton, m_newton, p_newton = ns.structure_solver('Newton', p0)
R_newton = r_newton[-1]
M_newton = m_newton[-1]

r_tov, m_tov, p_tov = ns.structure_solver('TOV', p0)
R_tov = r_tov[-1]
M_tov = m_tov[-1]

print('=========================================================')
print('Pure Neutron Star in the full relativistic  case.')
print('Central pressure : ', p0*conversion_dict['geom']['pressure']['cgs'], " dyne/cm^2")
print('---------------------------------------------------------')
print(' M_newton = ', M_newton, ' R_newton = ', R_newton,
      '\n', ' M_tov = ', M_tov, ' R_tov = ', R_tov)
print('=========================================================')

# create a figure of m(r) and P(r)

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.text(0.5, 1.07, "P(r) & m(r) of a pure neutron star in the full relativistic case",
         horizontalalignment='center',
         fontsize=12,
         transform = ax.transAxes)

ax.plot(r_newton, p_newton, color="blue", linestyle="-", linewidth=1, label = 'P Newton')
ax.plot(r_tov, p_tov, color="black", linestyle="-", linewidth=2,  label = 'P TOV')
ax.set_xlabel('r [km]',fontsize=10)
ax.set_ylabel(r'P [$dyne/cm^2$]', fontsize=10)
ax.minorticks_on()

ax2 = ax.twinx()
ax2.plot(r_newton, m_newton,color="blue", linestyle=":", label = 'm Newton')
ax2.plot(r_tov, m_tov, color="black", linestyle="-.", label = 'm TOV')
ax2.plot(R_newton, M_newton, marker = 'o', linestyle="", color='green', label='NS Newton mass')
ax2.plot(R_tov, M_tov, marker = 'o', linestyle="", color='red', label='NS TOV mass')
ax2.set_ylabel(r"m [$M_{\odot}$]",fontsize=10)
ax2.minorticks_on()

fig.legend(loc="upper center", bbox_to_anchor=(0.5,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'relns_mp-vs-r.pdf',
            format='pdf',
            dpi=1000)

#  Build a sequence of 500 stars. 

p0_1 = 3.5e32*conversion_dict['cgs']['pressure']['geom']
p0_2 =  3.5e38*conversion_dict['cgs']['pressure']['geom']

p_range = np.linspace(p0_1, p0_2, 200)

R_star_tov, M_star_tov = ns.mass_vs_radius('TOV', p_range)
R_star_newton, M_star_newton = ns.mass_vs_radius('Newton', p_range)

# plot Mass-Radius

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass-Radius of a pure neutron star in the full relativistic case")
ax.plot(R_star_newton, M_star_newton, color="blue", linestyle=":", linewidth=1, label = 'Newton')
ax.plot(R_star_tov, M_star_tov, color="black", linestyle=":", linewidth=2,  label = 'TOV')
ax.set_xlabel('R [km]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)
ax.minorticks_on()

fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'relns_mass-vs-radius.pdf',
            format='pdf',
            dpi=1000)

# plot Mass/Radius - central pressure 

fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.title("Mass/Radius vs Central Pressure of a pure NS in the full relativistic case")
ax.set_xscale('log')
ax.minorticks_on()
ax.plot(p_range*conversion_dict['geom']['pressure']['cgs'], M_star_newton, color="blue", linestyle="-.", linewidth=1, label = 'M-Newton')
ax.plot(p_range*conversion_dict['geom']['pressure']['cgs'], M_star_tov, color="black", linestyle="-.", linewidth=2,  label = 'M-TOV')
ax.set_xlabel('p0 [$dyne/cm^2$]',fontsize=14)
ax.set_ylabel(r"M [$M_{\odot}$]", fontsize=14)

ax2 = ax.twinx()
ax2.set_xscale('log')
ax2.minorticks_on()
ax2.plot(p_range*conversion_dict['geom']['pressure']['cgs'], R_star_newton, color="blue", linestyle=":", linewidth=1, label = 'R-Newton')
ax2.plot(p_range*conversion_dict['geom']['pressure']['cgs'], R_star_tov, color="black", linestyle=":", linewidth=2,  label = 'R-TOV')
ax2.set_ylabel(r"R [$km$]",fontsize=14)

fig.legend(loc='center right', bbox_to_anchor=(1,0.5), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"
fig.savefig(results_dir+'relns_mr-vs-p0.pdf',
            format='pdf',
            dpi=1000)