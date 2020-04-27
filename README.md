## Laser Plasma Interaction Cheat-Sheet

A simple script to calculate basic parameters for relativistic laser plasma interactions (for 
personal use mainly). 

### Motivation

When describing a laser, LPA people often define it by _energy on target_, or _output power_, 
with FWHM durations and beam size (can mean radius, waist, diameter etc). LaserPlasmaICS allows 
to define laser with different input. In the following example same laser is defined with the help 
of different parameters:

```python
from LPICS.ics import LaserPlasmaICS

lpi = LaserPlasmaICS(Energy=1.2, R_fwhm=16, tau_fwhm=27)
print(f"Power is {lpi.prm['Power']*1e-12:0.4g} TW" )
print(f"Duration field sqrt(2)*RMS is {lpi.prm['tau']:0.4g} fs" )

lpi = LaserPlasmaICS(Power=41.75e12, R_fwhm=16, tau=22.93)
print(f"\nField amplitude (normalized) is {lpi.prm['a0']:0.4g}" )
print(f"Beam waist is {lpi.prm['w0']:0.4g} um" )

lpi = LaserPlasmaICS(a0=2.595, w0=13.59, tau=22.93)
print(f"\nLaser energy is {lpi.prm['Energy']:0.4g} J")
```
The `lpi.prm` dictionary containes all necessary conversions. 

Besides conversions, for a given laser there are methods to calculate useful plasma densities 
corresponding to _matching_ and _critical_ conditions.
```python
# Densities for transverse (blowout) and longitudonal (linear wakefield)
print(f"{lpi.density_match('WLu'):0.4g}, {lpi.density_match('longitudinal'):0.4g}")

# Density at which laser power becomes critical
print(f"{lpi.density_match('critPower'):0.4g}")
```

### Installation

Can be installed by cloning the source 
```bash
git clone https://github.com/hightower8083/LaserPlasmaICS.git
cd LaserPlasmaICS
python setup.py install
```
or via PiPy
```bash
pip install git+https://github.com/hightower8083/LaserPlasmaICS.git
```

### Contributions

Feel free to propose your favorite formulas or fixes
