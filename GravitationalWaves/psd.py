"""Module contains functions to compute various power spectral densities for sensitivity curves."""

import numpy as np
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import splev, splrep
from importlib import resources

__all__ = ['load_response_function', 'approximate_response_function', 'power_spectral_density',
           'lisa_psd', 'get_confusion_noise', 'get_confusion_noise_robson19']

def load_response_function(f, fstar=19.09e-3):
    """Load in LISA response function from file and interpolate values for a range of frequencies.
    Parameters:
        f (float/array): Frequencies at which to evaluate the sensitivity curve
        fstar (float, optional): f* from Robson+19. Defaults to 19.09e-3 (19.09 mHz).
    Returns:
        R (function): LISA response function at each frequency.
    """    
    # try to load the values for interpolating R
    try:
        with resources.path(package="legwork", resource="R.npy") as path:
            f_R, R = np.load(path)
    except FileExistsError:  # pragma: no cover
        print("WARNING: Can't find response function file, using approximation instead")
        return approximate_response_function(f, fstar)

    # interpolate the R values in the file
    R_data = splrep(f_R * fstar, R, s=0)

    # use interpolated curve to get R values for supplied f values
    R = splev(f, R_data, der=0)
    return R

def approximate_response_function(f, fstar):
    """Approximate LISA response function
    Parameters:
        f (float/array): Frequencies at which to evaluate the sensitivity curve
        fstar ([type]): [description]
    """    
    return (3 / 10) / (1 + 0.6 * (f / fstar)**2)


def lisa_psd(f, t_obs=4 * u.yr, L=2.5e9 * u.m, approximate_R=False, confusion_noise="robson19"):
    """Calculates the effective LISA power spectral density sensitivity curve
    Parameters:
        f (float/array): Frequencies at which to evaluate the sensitivity curve
        t_obs ([type], optional): [description]. Defaults to 4*u.yr.
        L ([type], optional): [description]. Defaults to 2.5e9*u.m.
        approximate_R (bool, optional): [description]. Defaults to False.
        confusion_noise (str, optional): [description]. Defaults to "robson19".
    Returns:
        [type]: [description]
    """    
    # convert frequency from Hz to float for calculations
    f = f.to(u.Hz).value

    # minimum and maximum frequencies in Hz based on the R file from Robson+19
    MIN_F = 1e-7
    MAX_F = 2e0

    # overwrite frequencies that outside the range to prevent error
    f = np.where(np.logical_and(f >= MIN_F, f <= MAX_F), f, 1e-8)

    # single link optical metrology noise (Robson+ Eq. 10)
    def Poms(f):
        return (1.5e-11)**2 * (1 + (2e-3 / f)**4)

    # single test mass acceleration noise (Robson+ Eq. 11)
    def Pacc(f):
        return (3e-15)**2 * (1 + (0.4e-3 / f)**2) * (1 + (f / (8e-3))**4)

    # calculate response function (either exactly or with approximation)
    fstar = (const.c / (2 * np.pi * L)).to(u.Hz).value
    if approximate_R:
        R = approximate_response_function(f, fstar)
    else:
        R = load_response_function(f, fstar)

    # get confusion noise
    if isinstance(confusion_noise, str) or confusion_noise is None:
        cn = get_confusion_noise(f=f * u.Hz, t_obs=t_obs, model=confusion_noise).value
    else:
        cn = confusion_noise(f, t_obs).value
        
    L = L.to(u.m).value

    # calculate sensitivity curve
    psd = (1 / (L**2) * (Poms(f) + 4 * Pacc(f) / (2 * np.pi * f)**4)) / R + cn
    psd = np.where(np.logical_and(f >= MIN_F, f <= MAX_F), psd, np.inf)
    
    return psd / u.Hz



def power_spectral_density(f, instrument="LISA", custom_psd=None, t_obs="auto", L="auto",
                           approximate_R=False, confusion_noise="auto"):
    """Calculates the effective power spectral density for all instruments.
    Parameters:
        f (float/array): Frequencies at which to evaluate the sensitivity curve
        instrument (str, optional): [description]. Defaults to "LISA".
        custom_psd ([type], optional): [description]. Defaults to None.
        t_obs (str, optional): [description]. Defaults to "auto".
        L (str, optional): [description]. Defaults to "auto".
        approximate_R (bool, optional): [description]. Defaults to False.
        confusion_noise (str, optional): [description]. Defaults to "auto".
    Returns:
        [type]: [description]
    """    
    if instrument == "LISA":
        # update any auto values to be instrument specific
        L = 2.5e9 * u.m if L == "auto" else L
        confusion_noise = "robson19" if confusion_noise == "auto" else confusion_noise
        t_obs = 4 * u.yr if t_obs == "auto" else t_obs

        # calculate psd
        psd = lisa_psd(f=f, L=L, t_obs=t_obs, approximate_R=approximate_R, confusion_noise=confusion_noise)
    
    else:
        raise ValueError("instrument: `{}` not recognised".format(instrument))

    return psd


def get_confusion_noise_robson19(f, t_obs=4 * u.yr):
    """Calculate the confusion noise using the model from Robson+19 Eq. 14 and Table 1
    Parameters:
        f (float/array): Frequencies at which to evaluate the sensitivity curve
        t_obs ([type], optional): [description]. Defaults to 4*u.yr.
    Returns:
        confusion_noise (float): [description]
    """    
    # erase the units for speed
    f = f.to(u.Hz).value

    # parameters from Robson+ 2019 Table 1
    lengths = np.array([0.5, 1.0, 2.0, 4.0]) * u.yr
    alpha = [0.133, 0.171, 0.165, 0.138]
    beta = [243.0, 292.0, 299.0, -221.0]
    kappa = [482.0, 1020.0, 611.0, 521.0]
    gamma = [917.0, 1680.0, 1340.0, 1680.0]
    fk = [2.58e-3, 2.15e-3, 1.73e-3, 1.13e-3]

    # find index of the closest length to inputted observation time
    ind = np.abs(t_obs - lengths).argmin()

    confusion_noise = 9e-45 * f**(-7 / 3.) * np.exp(-f**(alpha[ind]) + beta[ind]
                                                    * f * np.sin(kappa[ind] * f))\
        * (1 + np.tanh(gamma[ind] * (fk[ind] - f))) * u.Hz**(-1)
    
    return confusion_noise


def get_confusion_noise(f, model, t_obs="auto"):
    """Calculate the confusion noise for a particular model
    Parameters:
        f (float/array): Frequencies at which to evaluate the sensitivity curve
        model (_type_): _description_
        t_obs (str, optional): _description_. Defaults to "auto".
    Returns:
        get_confusion_noise_robson19 (float): confusion noise
    """
    if model == "robson19":
        t_obs = 4 * u.yr if t_obs == "auto" else t_obs
        return get_confusion_noise_robson19(f=f, t_obs=t_obs)
    elif model is None:
        f = f.to(u.Hz).value
        return np.zeros_like(f) * u.Hz**(-1) if isinstance(f, (list, np.ndarray)) else 0.0 * u.Hz**(-1)
    else:
        raise ValueError("confusion noise model: `{}` not recognised".format(model))