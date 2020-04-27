# LaserPlasmaICS
Laser Plasma Interaction Cheat-Sheet

A simple script to calculate basic parameters for relativistic laser plasma interactions. 

```
from LPICS.ics import LaserPlasmaICS

lpi = LaserPlasmaICS(Energy=1.2, R_fwhm=16, tau_fwhm=27)

print(f"a0 = {lpi.prm['a0']:.2g}")
print(f"Power = {lpi.prm['Power'] * 1e-12:.2g} TW")
print(f"n_p (matched) = {lpi.density_match('WLu'):.2g} cm^-3")
print(f"n_p/n_c = {lpi.density_match('WLu')/lpi.density_match('crit'):.2g}")
```

For personal use mainly.