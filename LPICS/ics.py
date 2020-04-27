from scipy.constants import e, m_e, c, pi
from scipy.constants import epsilon_0 as eps0
import numpy as np

# Classic electron radius (meters)
r_e = e**2 / (4*pi * eps0 * m_e*c**2)

# Relativistic power unit (Watts)
P_ru = m_e * c**3 / r_e

# a0 for lambda=1um and I=1e18 W/cm^2
coef_I2a0 = (2/pi/P_ru * 1e10)**.5

# covert FWHM of `E^2` to RMS*sqrt(2) of `E`
coef_fwhm = (2*np.log(2))**-0.5

class LaserPlasmaICS:
    """
    Laser Plasma Interactive Cheat-Sheet

    After being initialized, has various conversions of input data
    in the dictionary `prm`

    Method:
    ---------
    density_match: finds values of plasma density for different
                   matching conditions (see documentation)
    """

    def __init__(l, **prm):
        """
        Initialize with given laser and plasma characteristics

        Parameter
        ---------
        lam0: micron
            Laser central wavelength. Default is set to 0.8um

        tau : fs
            laser duration defined as $E=E_0 \exp(-t^2/\tau^2)$.
            Alternatively, can be defined as full width at half maximum
            of power `tau_fwhm`

        w0: micron
            laser size at focus (waist) defined as $E=E_0 \exp(-r^2/w_0^2)$.
            Alternatively, can be defined as full width at half maximum
            of intensity `R_fwhm

        a0: unitless
            normalized vector-potential $a_0 = e E_max / (m_e c \omega_0)$
            Alternatively, can be defined with `Intensity`, `Power`, `Energy`
        """

        l.prm = prm

        if 'lam0' not in l.prm:
            l.prm['lam0'] = 0.8

        if 'tau_fwhm' in l.prm:
            l.prm['tau'] = l.prm['tau_fwhm'] * coef_fwhm
        else:
            l.prm['tau_fwhm'] = l.prm['tau'] / coef_fwhm

        if 'R_fwhm' in l.prm:
            l.prm['w0'] = l.prm['R_fwhm'] * coef_fwhm
        else:
            l.prm['R_fwhm'] = l.prm['w0'] / coef_fwhm

        if 'a0' in l.prm:
            l.prm['Intensity'] = 1e18 * (l.prm['a0']/l.prm['lam0']/coef_I2a0)**2
            l.prm['Power'] = l._Power_from_Intens()
            l.prm['Energy'] = l._Energy_from_Power()
        elif 'Intensity' in l.prm:
            l.prm['a0'] = l._a0_from_Intens()
            l.prm['Power'] = l._Power_from_Intens()
            l.prm['Energy'] = l._Energy_from_Power()
        elif 'Power' in l.prm:
            l.prm['Energy'] = l._Energy_from_Power()
            if 'w0' in l.prm:
                l.prm['Intensity'] = l._Intens_from_Power()
                l.prm['a0'] = l._a0_from_Intens()
        elif 'Energy' in l.prm:
            if 'tau' in l.prm:
                l.prm['Power'] = (2e30/pi)**.5 * l.prm['Energy']/l.prm['tau']
                if 'w0' in l.prm:
                    l.prm['Intensity'] = l._Intens_from_Power()
                    l.prm['a0'] = l._a0_from_Intens()

    def density_match(l, name):
        """
        Get density of electron plasma corresponding to
        one of avaliable matchings

        Parameter
        ---------
        name: string
            `WLu`         : transverse matching from [W Lu PRSTAB 2007]
            `crit`        : critical density for laser wavelength
            `longitudinal`: density for which laser duration is resonant
                            with a linear plasma wave [Gorbunov JETP 1987]
            `critPower`   : density for which laser power supports relativistic
                            self-focusing [G.-Z. Sun Phys. Fluids 1987]
        """
        if name is 'WLu':
            return 1e6 * l.prm['a0'] / pi  / r_e / l.prm['w0']**2
        if name is 'crit':
            return 1e6 * pi / r_e / l.prm['lam0']**2
        if name is 'longitudinal':
            return 1e6 / pi / r_e / (c*l.prm['tau']*1e-9)**2
        if name is 'critPower':
            n_pe = l.density_match('crit')
            return n_pe * (2*P_ru) / l.prm['Power']

    def _a0_from_Intens(l):
        if 'lam0' in l.prm:
            return 1e-9 * coef_I2a0*l.prm['lam0']*l.prm['Intensity']**.5
        else:
            print('Need laser wavelength')
            return 0.0

    def _Intens_from_Power(l):
        if 'w0' in l.prm:
            return 2e8/pi * l.prm['Power'] / l.prm['w0']**2
        else:
            print('Need laser waist')
            return 0.0

    def _Power_from_Intens(l):
        if 'w0' in l.prm:
            return pi/2e8 * l.prm['Intensity'] * l.prm['w0']**2
        else:
            print('Need laser waist')
            return 0.0

    def _Energy_from_Power(l):
        if 'tau' in l.prm:
            return (pi/2e30)**0.5 * l.prm['Power'] * l.prm['tau']
        else:
            print('Need laser duration')
            return 0.0