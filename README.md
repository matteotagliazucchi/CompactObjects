# compactobjects

A python package for studying astrophysical compact objects.

## Description

For the moment the package solves Tolman-Oppenheimer-Volkoff or Newton structures equations in order to determine the Mass and the Radius of astrophysical compact objects, such as neutron stars or white dwarfs. In order to solve such equations one should provide a proper equation of state (eos) and a value for the central pressure of the star.
The package works with two different type eos:
1. simple monotropic function in the form $$P = K\cdot\varepsilon^\gamma$$, managed by the class `PolytropicEos` which allows to return the energy density $$\varepsilon$$ in function of the pressure $$P$$. This type of eos is useful to study the non relativistic limit of a Fermi gas of a pure neutron stars or to study both the relativistic and the non relativistic limit of a pure white dwarf.
2. implicit form of the equations $$P = P(x)$$ and $$\varepsilon = \varepsilon(x)$$ where $$x$$ is the ratio between the Fermi momentum of the gas considered and $$m\cdot c$$, where $$m$$ is the mass of the particles of the gas. This case is managed by the class `ImplicitEos`, which generates a table of $$P$$ and $$\varepsilon$$ for a proper $$x$$ range. This class should be used for a pure neutron star in the full relativistic case, but also for the full relativistc case of a pure white dwarf. We can use this class also for studying neutron stars with the presence of protons and electrons.

## Prerequisites

1. numpy
2. scipy
3. matplotlib
4. astropy

## Installation

Place in your favorite folder.
Clone the repository in it

1. `git clone https://github.com/matteotagliazucchi/CompactObjects.git`

Then:

2. `cd CompactObjects`

If you want you can create a virtual environment (`virtual_env`) where install the package:

3. `python3 -m venv virtual_env`

where `virtual_env` can be replaced by `/path/to/directory/virtual_env` in the case you want the virtual envirnoment in the directory `/path/to/directory/`.
If you created a virtual environment, activate it:

4. `source virtual_env/bin/activate` or `source /path/to/directory/virtual_env/bin/activate`

Then proceed installing the package (the package will be installed only in the virtual environment, if you created one and you activated it). 
Installation steps are:

5. `python -m pip install -r ./requirements.txt`

to install the prerequisites and then

6. `python setup.py install`

or for installing in development mode:

7. `python setup.py develop`

If you want to install the package in your base python environemnt, ignore steps 3 and 4.

## Uninstallation

You can uninstall the package with the command:

1. `pip uninstall compactobjects`

once you have gone in the folder of the package `CompactObjects`
You can also eliminate the dependences contained in `requirements.txt` with the same command. However before uninstalling, you should ensure that the packages are not dependencies for other existing packages.

## Usage
See the folder 'tests' to see examples on how to use the package.

The basic steps are:
1. define an equation of state using `PolytropicEos` or `ImplicitEos`
2. define a `CompactStar` object by passing the eos instance to the constructor of the class
3. define a valure for the pressure at the center of the star. This value must be converted to geometrized units via the dictionary `conversion_dict['si']['pressure']['geom']` where the first 'si' can be replaced by 'cgs' i you want to provide the pressure in cgs units.
4. use the method `CompactStar.structuresolver()` which returns the mass $$m$$ and the pressure $$P$$ of the star in function of the distance from the center of the star $$r$$. The last items of these three returned array are the pressure at the surface of the star (~0) and thus the radius and the mass of the star. By default radius are given in km, pressures in dyne/cm^2 and mass in solar masses.
Then you can plot or do whatever you want.


## How to run tests

Follow installation steps. Eventually, remember to activate the virtual environment where you installed the package. 
Then, go inside the folder `CompactObjects/tests`. If you want to run a test on pure neutron stars, enter in the folder `PureNS`. Then open the terminal and run:

- `python nonrel_purens.py`

or 

- `python rel_purens.py`

Some results will be printed on the terminal, while plots will be available in the folder `NSOutput`

## References

For the physics used in this package see 
1. I. Sagert, M. Hempel, C. Greiner, J. Schaffner-Bielich; "*Compact Stars for Undergraduates*"; https://arxiv.org/abs/astro-ph/0506417v1
2. R.R. Silbar, S. Reddy; "*Neutron stars for undergraduates*"; https://doi.org/10.1119/1.1703544
