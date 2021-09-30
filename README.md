# compactobjects

A python package for studying astrophysical compact objects.

## Description

For the moment the package solves Tolman-Oppenheimer-Volkoff or Newton structures equations in order to determine the Mass and the Radius of astrophysical compact objects, such as neutron stars or white dwarfs. In order to solve such equations one should provide a proper equation of state (eos) and a value for the central pressure of the star.
The package works with different type of eos:
1. simple polytropic function in the form <img src="https://latex.codecogs.com/svg.image?P&space;=&space;K\cdot\varepsilon^\gamma" title="P = K\cdot\varepsilon^\gamma" />, managed by the class `PressureEdenPolytropic`. This type of eos is useful to study the non relativistic limit of a Fermi gas of a pure neutron stars or to study both the relativistic and the non relativistic limit of a pure white dwarf.
2. implicit form of the equations <img src="https://latex.codecogs.com/svg.image?P&space;=&space;P(x)" title="P = P(x)" /> and <img src="https://latex.codecogs.com/svg.image?\varepsilon&space;=&space;\varepsilon(x)" title="\varepsilon = \varepsilon(x)" /> where x is the ratio between the Fermi momentum of the gas considered and <img src="https://latex.codecogs.com/svg.image?m\cdot&space;c" title="m\cdot c" />, where m is the mass of the particles of the gas. This case is managed by the class `ImplicitEos`. This class should be used for a pure neutron star in the full relativistic case, but also for the full relativistc case of a pure white dwarf. We can use this class also for studying neutron stars with the presence of protons and electrons.
3. simple polytropic function in the form <img src="https://latex.codecogs.com/svg.image?P&space;=&space;K\cdot\rho^\gamma" title="P = K\cdot\rho^\gamma" />, managed by the class `PressureDensityPolytropic`. 
4. piecewise polytropic eos, where we have a simple polytropic functions in each of the user defined range. This kind of eos are implemented in the class `PressureDensityPiecewise`

The last two classes are used to study various piecewise polytropic eos's, classified by Read et al in [3].

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

See the folder 'tests' to see examples on how to use the package. The various tests script show how to use all the classe:
1. `nonrel_purens.py`, `nonrel_wd.py` and `rel_wd.py` use the eos class `PressureEdenPolytropic`
2. `rel_purens.py` uses `ImplicitEos`
3. `read_ns.py` basically uses `PressureDensityPiecewise`.

All these scripts produces:
- the profile m(r) and p(r) of the studied compact object in a particular value of the central pressure
- the mass-vs-radius profile of the compact object
- the mass-vs-central pressure and the radius-vs-centralpressure profile of the compact object

Sometimes scripts produce both the classical curve (solving Newton structure equation) and those with relativistic effects (TOV equations). When it's not specified, they solve only TOV equations.


## How to run tests

Follow installation steps. Eventually, remember to activate the virtual environment where you installed the package. 
Then, go inside the folder `CompactObjects/tests`. If you want to run a test on pure neutron stars, enter in the folder `PureNS`. Then open the terminal and run:

- `python nonrel_purens.py`

or 

- `python rel_purens.py`

Some results will be printed on the terminal, while plots will be available in the folder `NSOutput`.
A similar procedure holds also for the other tests.

## Acknowledgements

I would like to thank Janko Petković and Lorenzo Squadrani for helping me in the comprehension of virtual environments and of `setup.py` files.

## References

For the physics used in this package see 
1. I. Sagert, M. Hempel, C. Greiner, J. Schaffner-Bielich; "*Compact Stars for Undergraduates*"; https://arxiv.org/abs/astro-ph/0506417v1
2. R.R. Silbar, S. Reddy; "*Neutron stars for undergraduates*"; https://doi.org/10.1119/1.1703544
3. J.S. Read, B.D. Lackey, B.J. Owen, J.L. Friedman; "*Constraints on a phenomenologically parametrized neutron-star equation of state*"; https://doi.org/10.1103/PhysRevD.79.124032
