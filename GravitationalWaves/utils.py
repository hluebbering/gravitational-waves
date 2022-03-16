from scipy.special import jv
from astropy import constants as c
from astropy import units as u
import numpy as np

__all__ = ['chirp_mass', 'peters_g', 'peters_f', 'get_a_from_f_orb', 'get_f_orb_from_a', 'get_a_from_ecc',
           'beta', 'c_0', 'fn_dot', 'ensure_array', 'D_plus_squared',
           'D_cross_squared', 'D_plus_D_cross', 'F_plus_squared', 'F_cross_squared']


def get_a_from_f_orb(f_orb, m_1, m_2):
    """Converts orbital frequency to semi-major axis
    Using Kepler's third law, convert orbital frequency to semi-major axis.
    Inverse of :func:`legwork.utils.get_f_orb_from_a`.
    Parameters
    ----------
    f_orb : `float/array`
        Orbital frequency
    m_1 : `float/array`
        Primary mass
    m_2 : `float/array`
        Secondary mass
    Returns
    -------
    a : `float/array`
        Semi-major axis
    """
    a = (c.G * (m_1 + m_2) / (2 * np.pi * f_orb)**2)**(1/3)

    # simplify units if present
    if isinstance(a, u.quantity.Quantity):
        a = a.to(u.AU)

    return a


def get_f_orb_from_a(a, m_1, m_2):
    """Converts semi-major axis to orbital frequency
    Using Kepler's third law, convert semi-major axis to orbital frequency.
    Inverse of :func:`legwork.utils.get_a_from_f_orb`.
    Parameters
    ----------
    a : `float/array`
        Semi-major axis
    m_1 : `float/array`
        Primary mass
    m_2 : `float/array`
        Secondary mass
    Returns
    -------
    f_orb : `float/array`
        Orbital frequency
    """
    f_orb = ((c.G * (m_1 + m_2) / a**3))**(0.5) / (2 * np.pi)

    # simplify units if present
    if isinstance(f_orb, u.quantity.Quantity):
        f_orb = f_orb.to(u.Hz)

    return f_orb


def ensure_array(*args):
    """Convert arguments to numpy arrays
    Convert arguments based on the following rules
        - Ignore any None values
        - Convert any lists to numpy arrays
        - Wrap any other types in lists and convert to numpy arrays
    Parameters
    ----------
    args : `any`
        Supply any number of arguments of any type
    Returns
    -------
    array_args : `any`
        Args converted to numpy arrays
    any_not_arrays : `bool`
        Whether any arg is not a list or None or a numpy array
    """
    array_args = [None for i in range(len(args))]
    any_not_arrays = False
    for i in range(len(array_args)):
        exists = args[i] is not None
        has_units = isinstance(args[i], u.quantity.Quantity)
        if exists and has_units:
            if not isinstance(args[i].value, np.ndarray):
                any_not_arrays = True
                array_args[i] = np.asarray([args[i].value]) * args[i].unit
            else:
                array_args[i] = args[i]
        elif exists and not has_units:
            if not isinstance(args[i], np.ndarray):
                if not isinstance(args[i], list):
                    any_not_arrays = True
                    array_args[i] = np.asarray([args[i]])
                else:
                    array_args[i] = np.asarray(args[i])
            else:
                array_args[i] = args[i]
        else:
            array_args[i] = args[i]
    return array_args, any_not_arrays


def chirp_mass(m_1, m_2):
    """ 
    Computes chirp mass of binaries
    Parameters: 
        m_1 : Primary mass (`float/array`)
        m_2 : Secondary mass (`float/array`)
    Returns:
        m_c : Chirp mass (`float/array`)
    """
    m_c = (m_1 * m_2)**(3/5) / (m_1 + m_2)**(1/5)

    # simplify units if present
    if isinstance(m_c, u.quantity.Quantity):
        m_c = m_c.to(u.Msun)

    return m_c


def peters_g(n, e):
    """Compute g(n, e) from Peters and Mathews (1963) Eq.20
    This function gives the relative power of gravitational radiation at the nth harmonic
    Parameters
    ----------
    n : `int/array`
        Harmonic(s) of interest
    e : `float/array`
        Eccentricity
    Returns
    -------
    g : `array`
        g(n, e) from Peters and Mathews (1963) Eq. 20
    """

    bracket_1 = jv(n-2, n*e) - 2*e*jv(n-1, n*e) + 2/n*jv(n, n*e) + 2*e*jv(n+1, n*e) - jv(n+2, n*e)
    bracket_2 = jv(n-2, n*e) - 2*jv(n, n*e) + jv(n+2, n*e)
    bracket_3 = jv(n, n*e)

    g = n**4/32 * (bracket_1**2 + (1 - e**2) * bracket_2**2 + 4 / (3 * n**3) * bracket_3**2)

    return g


def peters_f(e):
    """f(e) from Peters and Mathews (1963) Eq.17
    This function gives the integrated enhancement factor of gravitational radiation from an eccentric
    source compared to an equivalent circular source.
    Parameters
    ----------
    e : `float/array`
        Eccentricity
    Returns
    -------
    f : `float/array`
        Enhancement factor
    Notes
    -----
    Note that this function represents an infinite sum of g(n, e) - :meth:`legwork.utils.peters_g`
    """

    numerator = 1 + (73/24)*e**2 + (37/96)*e**4
    denominator = (1 - e**2)**(7/2)

    f = numerator / denominator

    return f


