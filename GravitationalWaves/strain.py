# Computes several types of gravitational wave strains

import astropy.constants as c
from GravitationalWaves import utils
import numpy as np

__all__ = ['h_0_n', 'h_c_n']

# Computes strain amplitude at n harmonics
def h_0_n(m_c, f_orb, ecc, n, dist, position=None, polarisation=None, inclination=None, interpolated_g=None):
    """
    Parameters
    m_c : `float/array` Chirp mass of each binary.
    f_orb : `float/array` Orbital frequency of each binary at each timestep.
    ecc : `float/array` Eccentricity of each binary at each timestep. 
    n : `int/array` Harmonic(s) at which to calculate the strain. 
    dist : `float/array` Distance to each binary. 
    
    position : `SkyCoord/array` Sky position of source. 
    polarisation : `float/array` GW polarisation angle of the source. 
    inclination : `float/array` Inclination of the source. 
    interpolated_g : `function` Computes g(n,e) from Peters (1964).
    Returns
    h_0 : `float/array` Strain amplitude. Shape is (x, y, z).
    """
    # convert to array if necessary
    arrayed_args, _ = utils.ensure_array(m_c, f_orb, ecc, n, dist)
    m_c, f_orb, ecc, n, dist = arrayed_args

    # if one timestep then extend dimensions
    if f_orb.ndim != 2:
        f_orb = f_orb[:, np.newaxis]
    if ecc.ndim != 2:
        ecc = ecc[:, np.newaxis]

    # extend mass and distance dimensions
    m_c = m_c[:, np.newaxis]
    dist = dist[:, np.newaxis]

    # work out strain for n independent part and broadcast to correct shape
    prefac = (2**(28/3) / 5)**(0.5) * c.G**(5/3) / c.c**4
    n_independent_part = prefac * m_c**(5/3) * (np.pi * f_orb)**(2/3) / dist

    # check whether to interpolate g(n, e)
    if interpolated_g is None:
        # extend harmonic and eccentricity dimensions to full (x, y, z)
        n = n[np.newaxis, np.newaxis, :]
        ecc = ecc[..., np.newaxis]
        n_dependent_part = utils.peters_g(n, ecc)**(1/2) / n
    else:
        # flatten array to work nicely interp2d
        g_vals = interpolated_g(n, ecc.flatten())

        # set negative values from cubic fit to 0.0
        g_vals[g_vals < 0.0] = 0.0

        # unsort the output array if there is more than one eccentricity
        if isinstance(ecc, (np.ndarray, list)) and len(ecc) > 1:
            g_vals = g_vals[np.argsort(ecc.flatten()).argsort()]

        # reshape output to proper dimensions
        g_vals = g_vals.reshape((*ecc.shape, len(n)))

        n_dependent_part = g_vals**(0.5) / n[np.newaxis, np.newaxis, :]

    h_0 = n_independent_part[..., np.newaxis] * n_dependent_part

    if position is not None or polarisation is not None or inclination is not None:
        raise ValueError("Specific source position, polarisation, and/or inclination not available.")

    return h_0.decompose().value

# Computes strain amplitude at n harmonics
def h_c_n(m_c, f_orb, ecc, n, dist, position=None, polarisation=None, inclination=None, interpolated_g=None):
    
    # convert to array if necessary
    arrayed_args, _ = utils.ensure_array(m_c, f_orb, ecc, n, dist)
    m_c, f_orb, ecc, n, dist = arrayed_args

    # if one timestep then extend dimensions
    if f_orb.ndim != 2:
        f_orb = f_orb[:, np.newaxis]
    if ecc.ndim != 2:
        ecc = ecc[:, np.newaxis]

    # extend mass and distance dimensions
    m_c = m_c[:, np.newaxis]
    dist = dist[:, np.newaxis]

    # work out strain for n independent part
    prefac = (2**(5/3) / (3 * np.pi**(4/3)))**(0.5) * c.G**(5/6) / c.c**(3/2)
    n_independent_part = prefac * m_c**(5/6) / dist * f_orb**(-1/6) / utils.peters_f(ecc)**(0.5)

    # check whether to interpolate g(n, e)
    if interpolated_g is None:
        # extend harmonic and eccentricity dimensions to full (x, y, z)
        n = n[np.newaxis, np.newaxis, :]
        ecc = ecc[..., np.newaxis]
        n_dependent_part = (utils.peters_g(n, ecc) / n)**(1/2)
    else:
        # flatten array to work nicely interp2d
        g_vals = interpolated_g(n, ecc.flatten())

        # set negative values from cubic fit to 0.0
        g_vals[g_vals < 0.0] = 0.0

        # unsort the output array if there is more than one eccentricity
        if isinstance(ecc, (np.ndarray, list)) and len(ecc) > 1:
            g_vals = g_vals[np.argsort(ecc.flatten()).argsort()]

        # reshape output to proper dimensions
        g_vals = g_vals.reshape((*ecc.shape, len(n)))

        n_dependent_part = (g_vals / n[np.newaxis, np.newaxis, :])**(0.5)

    h_c = n_independent_part[..., np.newaxis] * n_dependent_part

    if position is not None or polarisation is not None or inclination is not None:
        raise ValueError("Specific source position, polarisation, and/or inclination not available.")

    return h_c.decompose().value