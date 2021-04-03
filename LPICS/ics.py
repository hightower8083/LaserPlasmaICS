from scipy.constants import e, m_e, c, pi
from scipy.constants import epsilon_0 as eps0
from scipy.constants import k as k_B
import numpy as np

try:
    from mendeleev import element as table_element
    mendeleev_avail = True
except Exception:
    mendeleev_avail = False


# Classic electron radius (meters)
r_e = e**2 / (4*pi * eps0 * m_e*c**2)

# Relativistic power unit (Watts)
P_ru = m_e * c**3 / r_e

# a0 for lambda=1um and I=1e18 W/cm^2
coef_I2a0 = (2/pi/P_ru * 1e10)**.5

# covert FWHM of `E^2` to RMS*sqrt(2) of `E`
coef_fwhm = (2*np.log(2))**-0.5

mc2_MeV = m_e * c**2 / e * 1e-6

m_e_cgs = 9.1093837015e-28
c_cgs = 29979245800.0
e_cgs = 4.803204712570263e-10
coef_Eion2a0 = 1e-4 * (e * 1e7)**2 / ( 8*pi * m_e_cgs * c_cgs**2 * e_cgs**2)

class CheatSheet:
    """
    Laser Plasma Interaction Cheat-Sheet

    After being initialized, has various conversions of input data
    in the dictionary `prm`

    Method:
    ---------
    match_density: finds values of plasma density for different
                   matching conditions (see documentation)
    """

    def __init__(l, **prm):
        """
        Initialize with given laser and plasma characteristics

        Parameter
        ---------
        lam0: micron
            Laser central wavelength in microns. Default is set to 0.8 um

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

        n_pe: cm^-3
            electron plasma density

        T_e: keV
            electron plasma temperature
        """

        l.prm = prm

        if 'verbose' not in l.prm:
            l.verbose = False
        else:
            l.verbose = l.prm['verbose']

        # LASER
        if 'lam0' not in l.prm:
            l.prm['lam0'] = 0.8
            if l.verbose: print('assume lam0=0.8')

        l.prm['k0'] = 2 * np.pi / l.prm['lam0']
        l.prm['n_c'] = 1e6 * pi / r_e / l.prm['lam0']**2

        if 'pol' not in l.prm:
            l.prm['pol'] = 'linear'
            if l.verbose: print("assume pol='linear'")

        if 'tau' in l.prm:
            l.prm['ctau'] = l.prm['tau'] * c * 1e-9
            l.prm['tau_fwhm'] = l.prm['tau'] / coef_fwhm
            l.prm['ctau_fwhm'] = l.prm['tau_fwhm'] * c * 1e-9
        elif 'tau_fwhm' in l.prm:
            l.prm['ctau_fwhm'] = l.prm['tau_fwhm'] * c * 1e-9
            l.prm['tau'] = l.prm['tau_fwhm'] * coef_fwhm
            l.prm['ctau'] = l.prm['tau'] * c * 1e-9
        else:
            if l.verbose: print('no pulse duration')

        if 'w0' in l.prm:
            l.prm['R_fwhm'] = l.prm['w0'] / coef_fwhm
            l.prm['Rayleigh'] = np.pi * l.prm['w0']**2 / l.prm['lam0']
        elif 'R_fwhm' in l.prm:
            l.prm['w0'] = l.prm['R_fwhm'] * coef_fwhm
            l.prm['Rayleigh'] = np.pi * l.prm['w0']**2 / l.prm['lam0']
        else:
            if l.verbose: print('no pulse size')

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
            if l.verbose: print('no laser field')

        if 'a0' in l.prm:
            if l.prm['pol'] == 'linear':
                 l.prm['gamma_p'] = (1 + l.prm['a0']**2/2.)**.5
            elif l.prm['pol'] == 'circular':
                 l.prm['gamma_p'] = (1 + l.prm['a0']**2)**.5

        # PLASMA
        if 'T_e' in l.prm:
            l.prm['v_e'] = (e * l.prm['T_e'] / m_e)**.5 / c
            l.prm['T_e_K'] = e * l.prm['T_e'] / k_B

            if l.prm['v_e']>1.0:
                l.prm['v_e'] = 1.0
                if l.verbose: print('plasma is relativistic')
        else:
            l.prm['T_e'] = np.infty
            l.prm['T_e_K'] = np.infty
            l.prm['v_e'] = 1.0

        if 'n_pe' in l.prm:
            l.prm['k_pe0'] = (4*pi * r_e * l.prm['n_pe']*1e6)**.5
            l.prm['omega_pe0'] = l.prm['k_pe0'] * c
            l.prm['lam_De0'] = l.prm['v_e']/l.prm['k_pe0'] * 1e6

            if 'gamma_p' in l.prm:
                l.prm['n_pe_rel'] = l.prm['n_pe'] / l.prm['gamma_p']
                l.prm['k_pe'] = l.prm['k_pe0'] / l.prm['gamma_p']**.5
                l.prm['omega_pe'] = l.prm['k_pe'] * c
                l.prm['lam_De'] = l.prm['v_e']/l.prm['k_pe'] * 1e6
                l.prm['N_pe'] = l.prm['n_pe_rel']/l.prm['n_c']
            else:
                l.prm['k_pe'] = l.prm['k_pe0']
                l.prm['omega_pe'] = l.prm['omega_pe0']
                l.prm['lam_De'] = l.prm['lam_De0']
                l.prm['N_pe'] = l.prm['n_pe']/l.prm['n_c']

            if l.prm['N_pe'] < 1.0:
                l.prm['v_wake'] = (1 + l.prm['N_pe'])**-0.5

                if 'w0' in l.prm:
                    l.prm['v_wake'] -= l.prm['lam0']**2 / (2*pi*l.prm['w0'])**2
            else:
                l.prm['v_wake'] = 0.0

            l.prm['gamma_w'] = ( 1 - l.prm['v_wake']**2 )**-0.5

    def match_density(l, name):
        """
        Get density of electron plasma corresponding to
        one of avaliable matchings

        Parameter
        ---------
        name: string
            'WLu'         : transverse matching from [W Lu PRSTAB 2007]
            'transverse'  : same as `WLu` but with $a_0$ replaced by
                            $\gamma_p\sqrt{2}$
            'longitudinal': density for which laser duration is resonant
                            with a linear plasma wave [Gorbunov JETP 1987]
            'longitud_rel': same as `longitudinal` but with added Lorentz
                            factor $\gamma_p$
            'critPower'   : density for which laser power supports relativistic
                            self-focusing [G.-Z. Sun Phys. Fluids 1987]
        """
        if name == 'WLu':
            return 1e6 * l.prm['a0'] / pi  / r_e / l.prm['w0']**2
        if name == 'transverse':
            return 1e6 * 2**.5 * l.prm['gamma_p'] / pi  / r_e / l.prm['w0']**2
        if name == 'longitudinal':
            return 1e6 / pi / r_e / (c*l.prm['tau']*1e-9)**2
        if name == 'longitud_rel':
            return 1e6 * l.prm['gamma_p'] / pi / r_e / (c*l.prm['tau']*1e-9)**2
        if name == 'critPower':
            return l.prm['n_c'] * (2*P_ru) / l.prm['Power']

    def element_Zmax(l, name):
        """
        Get the maximal ion number which can be produced by the laser
        optical field ionization. The expression is taken from a simple
        static field model given in [S. Augst et al PRL 63 2212 (1989)].

        Parameter
        ---------
        name: string
            Name of the element which can be interpreted by medeleev library
        """
        if not mendeleev_avail:
            print('mendeleev library is not avaliable')
            return 0.0

        ionenergies = table_element(name).ionenergies
        Z_states = np.arange(1, len(ionenergies)+1)

        a_ioniz = np.array([ionenergies[iz]**2/iz for iz in Z_states])
        a_ioniz *= coef_Eion2a0*l.prm['lam0']

        return (a_ioniz<l.prm['a0']).sum()

    def lwfa_max_energy(l, name):
        """
        Get maximum electron energy achieved in LWFA process.

        Parameter
        ---------
        name: string
            'WLu'         : first part in Eq (6).a from [W Lu PRSTAB 2007]
            'WLu1'        : second part from Eq (6).a from [W Lu PRSTAB 2007]
            'Gordienko'   : Eq (1) from [Gordienko, Pukhov PoP 2005]
            'GordienkoLu' : Eq (13).a from [W Lu PRSTAB 2007]
        """
        if name == 'WLu':
            return 2/3 * mc2_MeV * l.prm['n_c']/ l.prm['n_pe'] * l.prm['a0']
        if name == 'WLu1':
            return mc2_MeV * (l.prm['n_c']/ l.prm['n_pe'])**(2/3) \
                    * (l.prm['Power'] / P_ru)**(1/3)
        if name == 'Gordienko':
            return 0.65 * mc2_MeV * (l.prm['ctau_fwhm']/l.prm['lam0']) \
                    * (l.prm['Power'] / P_ru)**0.5
        if name == 'GordienkoLu':
            return 0.16 * mc2_MeV * (l.prm['ctau_fwhm']/l.prm['w0']) \
                    * (l.prm['Power'] / P_ru)**(2/3) \
                    * (l.prm['n_c']/ l.prm['n_pe'])**(1/3)

    def lwfa_max_charge(l, name):
        """
        Get maximum electron beam charge achieved in LWFA process (in pC).

        Parameter
        ---------
        name: string
            'WLu'         : Eq (10).a from [W Lu PRSTAB 2007]
            'Gordienko'   : Eq (2) from [Gordienko, Pukhov PoP 2005]
        """
        if name == 'WLu':
            return 0.53 * e * 1e12 / (l.prm['k0'] * r_e * 1e6) \
                    * (l.prm['Power'] / P_ru)**(1/2)
        if name == 'Gordienko':
            return 1.8 * e * 1e12 / (l.prm['k0'] * r_e * 1e6) \
                    * (l.prm['Power'] / P_ru)**(1/2)

    def _a0_from_Intens(l):
        return 1e-9 * coef_I2a0*l.prm['lam0']*l.prm['Intensity']**.5

    def _Intens_from_Power(l):
        if 'w0' in l.prm:
            return 2e8/pi * l.prm['Power'] / l.prm['w0']**2
        else:
            if l.verbose:  print('Need laser waist, Intensity is set to 0')
            return 0.0

    def _Power_from_Intens(l):
        if 'w0' in l.prm:
            return pi/2e8 * l.prm['Intensity'] * l.prm['w0']**2
        else:
            if l.verbose:  print('Need laser waist, Intensity is set to 0')
            return 0.0

    def _Energy_from_Power(l):
        if 'tau' in l.prm:
            return (pi/2e30)**0.5 * l.prm['Power'] * l.prm['tau']
        else:
            if l.verbose: print('Need pulse duration, Energy is set to 0')
            return 0.0
