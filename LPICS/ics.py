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
    Laser Plasma Interaction Cheat-Sheet

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
            Laser central wavelength in microns. Default is set to 0.8

        tau : fs
            laser duration defined as $E=E_0 \exp(-t^2/\tau^2)$.
            Alternatively, can be defined as full width at half maximum
            of power `tau_fwhm`

        w0: micron
            laser size at focus (waist) defined as $E=E_0 \exp(-r^2/w_0^2)$.
            Alternatively, can be defined as full width at half maximum
            of intensity `R_fwhm`

        a0: unitless
            normalized vector-potential $a_0 = e E_max / (m_e c \omega_0)$
            Alternatively, can be defined with `Intensity` [W/cm^2],
            `Power` [W], `Energy` [J]

        pol: `linear` (default) or `circular`
            laser polarization. Accounted in the plasma Lorentz factor
            $\gamma_p$ being $\sqrt{1+a_0^2}$ for `circular`, and
            $\sqrt{1+a_0^2/2}$ for `linear` polarisations
        """

        l.prm = prm

        if 'verbose' not in l.prm:
            verbose = False
        else:
            verbose = l.prm['verbose']

        if 'lam0' not in l.prm:
            l.prm['lam0'] = 0.8
            if verbose: print('assume lam0=0.8')

        if 'pol' not in l.prm:
            l.prm['pol'] = 'linear'
            if verbose: print("assume pol='linear'")

        if 'tau_fwhm' in l.prm:
            l.prm['tau'] = l.prm['tau_fwhm'] * coef_fwhm
        elif 'tau' in l.prm:
            l.prm['tau_fwhm'] = l.prm['tau'] / coef_fwhm
        else:
            if verbose: print('no pulse duration')

        if 'R_fwhm' in l.prm:
            l.prm['w0'] = l.prm['R_fwhm'] * coef_fwhm
        elif 'w0' in l.prm:
            l.prm['R_fwhm'] = l.prm['w0'] / coef_fwhm
        else:
            if verbose: print('no pulse size')

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
        else:
            if verbose: print('no laser field')

        if 'a0' in l.prm:
            if l.prm['pol'] is 'linear':
                 l.prm['gamma_p'] = (1 + l.prm['a0']**2/2.)**.5
            elif l.prm['pol'] is 'circular':
                 l.prm['gamma_p'] = (1 + l.prm['a0']**2)**.5

    def density_match(l, name):
        """
        Get density of electron plasma corresponding to
        one of avaliable matchings

        Parameter
        ---------
        name: string
            `crit`        : critical density for laser wavelength
            `WLu`         : transverse matching from [W Lu PRSTAB 2007]
            `transverse`  : same as `WLu` but with $a_0$ replaced by
                            $\gamma_p\sqrt{2}$
            `longitudinal`: density for which laser duration is resonant
                            with a linear plasma wave [Gorbunov JETP 1987]
            `longitud_rel`: same as `longitudinal` but with added Lorentz
                            factor $\gamma_p$
            `critPower`   : density for which laser power supports relativistic
                            self-focusing [G.-Z. Sun Phys. Fluids 1987]
        """
        if name is 'crit':
            return 1e6 * pi / r_e / l.prm['lam0']**2
        if name is 'WLu':
            return 1e6 * l.prm['a0'] / pi  / r_e / l.prm['w0']**2
        if name is 'transverse':
            return 1e6 * 2**.5 * l.prm['gamma_p'] / pi  / r_e / l.prm['w0']**2
        if name is 'longitudinal':
            return 1e6 / pi / r_e / (c*l.prm['tau']*1e-9)**2
        if name is 'longitud_rel':
            return 1e6 * l.prm['gamma_p'] / pi / r_e / (c*l.prm['tau']*1e-9)**2
        if name is 'critPower':
            n_pe = l.density_match('crit')
            return n_pe * (2*P_ru) / l.prm['Power']

    def _a0_from_Intens(l):
        return 1e-9 * coef_I2a0*l.prm['lam0']*l.prm['Intensity']**.5

    def _Intens_from_Power(l):
        if 'w0' in l.prm:
            return 2e8/pi * l.prm['Power'] / l.prm['w0']**2
        else:
            if verbose:  print('Need laser waist, Intensity is set to 0')
            return 0.0

    def _Power_from_Intens(l):
        if 'w0' in l.prm:
            return pi/2e8 * l.prm['Intensity'] * l.prm['w0']**2
        else:
            if verbose:  print('Need laser waist, Intensity is set to 0')
            return 0.0

    def _Energy_from_Power(l):
        if 'tau' in l.prm:
            return (pi/2e30)**0.5 * l.prm['Power'] * l.prm['tau']
        else:
            if verbose: print('Need pulse duration, Energy is set to 0')
            return 0.0