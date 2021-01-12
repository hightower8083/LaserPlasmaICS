[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/hightower8083/LaserPlasmaICS.git/master?filepath=.%2Fbinder_example.ipynb)

## Laser Plasma Interaction Cheat-Sheet

A simple script to calculate basic parameters for relativistic laser plasma interactions. 

### Motivation

When describing a laser, LPA people often define it by _energy on target_, or _output power_, 
with FWHM durations and beam size (can mean radius, waist, diameter etc). CheatSheet is a tiny library, which allows 
to define laser with different inputs and provides convertions between them. In the following example same laser is defined with the help 
of different parameters:

```python
from LPICS.ics import CheatSheet as cs

lpi = cs(Energy=1.2, R_fwhm=16, tau_fwhm=27)
print(f"Power is {lpi.prm['Power']*1e-12:0.4g} TW" )
print(f"Duration field sqrt(2)*RMS is {lpi.prm['tau']:0.4g} fs" )

lpi = cs(Power=41.75e12, R_fwhm=16, tau=22.93)
print(f"\nField amplitude (normalized) is {lpi.prm['a0']:0.4g}" )
print(f"Beam waist is {lpi.prm['w0']:0.4g} um" )

lpi = cs(a0=2.595, w0=13.59, tau=22.93)
print(f"\nLaser energy is {lpi.prm['Energy']:0.4g} J")
```

```
Power is 41.75 TW
Duration field sqrt(2)*RMS is 22.93 fs

Field amplitude (normalized) is 2.595
Beam waist is 13.59 um

Laser energy is 1.2 J
```

The `lpi.prm` dictionary containes the necessary conversions, and can be used directly, e.g. `cs(Power=1e15, R_fwhm=50).prm['a0']`, `cs(lam0=0.8).prm['n_c']`. Besides conversions, for a given laser there are useful methods to calculate some _matching_ plasma densities:
-  `lpi.match_density('WLu')`
-  `lpi.match_density('longitudinal')`
-  `lpi.match_density('critPower')`

or to estimate the maximum ionization level produced by the laser:
-  `lpi.element_Zmax(name=element_name, e.g. 'Ar')`

### Try online

You can use the script [online via Binder](https://mybinder.org/v2/gh/hightower8083/LaserPlasmaICS.git/master?filepath=.%2Fbinder_example.ipynb)

### Installation

Can be installed by cloning the source 
```bash
git clone https://github.com/hightower8083/LaserPlasmaICS.git
cd CheatSheet
pip install .
```
or via PiPy
```bash
pip install git+https://github.com/hightower8083/LaserPlasmaICS.git
```

### Contributions

Feel free to propose your favorite formulas or fixes
