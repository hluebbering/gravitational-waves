"""
Module for computing several types of gravitational wave strains.
"""

import astropy.constants as c
from GravitationalWaves import utils
import numpy as np


__all__ = ['h_0_n', 'h_c_n']

def h_0_n(m_c, f_orb, ecc, n, dist, position=None, polarisation=None, inclination=None, interpolated_g=None):
    
    """
    Computes strain amplitude at n harmonics
    
    ### Parameters:
        m_c (`float/array`): Chirp mass of each binary.
        f_orb (`float/array`): Orbital frequency of each binary at each timestep.
        ecc (`float/array`): Eccentricity of each binary at each timestep. 
        n (`int/array`): Harmonic(s) at which to calculate the strain. 
        dist (`float/array`): Distance to each binary. 
        position (`SkyCoord/array`, optional): Sky position of source. Defaults to None.
        polarisation (`float/array`, optional): GW polarisation angle of the source. Defaults to None.
        inclination (`float/array`, optional): Inclination of the source. Defaults to None.
        interpolated_g (`function`, optional): Computes g(n,e) from Peters (1964). Defaults to None.
    
    ### Returns:
        h_0 (`float/array`): Strain amplitude. Shape is (x, y, z).
    """
    
    arrayed_args, _ = utils.ensure_array(m_c, f_orb, ecc, n, dist) # convert to array
    m_c, f_orb, ecc, n, dist = arrayed_args
    
    if f_orb.ndim != 2:  # extend dimensions
        f_orb = f_orb[:, np.newaxis]
    if ecc.ndim != 2:
        ecc = ecc[:, np.newaxis]

    m_c = m_c[:, np.newaxis] # extend mass and distance dimensions
    dist = dist[:, np.newaxis]

    # get strain for each n and convert to correct shape
    prefac = (2**(28/3) / 5)**(0.5) * c.G**(5/3) / c.c**4
    n_independent_part = prefac * m_c**(5/3) * (np.pi * f_orb)**(2/3) / dist

    if interpolated_g is None: # check whether to interpolate g(n, e)
        n = n[np.newaxis, np.newaxis, :] # extend harmonic, eccentricity dimensions
        ecc = ecc[..., np.newaxis]
        n_dependent_part = utils.peters_g(n, ecc)**(1/2) / n
    else:
        g_vals = interpolated_g(n, ecc.flatten()) # flatten array for interp2d
        g_vals[g_vals < 0.0] = 0.0 # set negative values from cubic fit to 0
        
        if isinstance(ecc, (np.ndarray, list)) and len(ecc) > 1: # unsort output array
            g_vals = g_vals[np.argsort(ecc.flatten()).argsort()]
            
        g_vals = g_vals.reshape((*ecc.shape, len(n))) # reshape output dimensions
        n_dependent_part = g_vals**(0.5) / n[np.newaxis, np.newaxis, :]

    h_0 = n_independent_part[..., np.newaxis] * n_dependent_part
    if position is not None or polarisation is not None or inclination is not None:
        raise ValueError("Specific source position, polarisation, and/or inclination not available.")

    return h_0.decompose().value


def h_c_n(m_c, f_orb, ecc, n, dist, position=None, polarisation=None, inclination=None, interpolated_g=None):
    
    """
    Computes strain amplitude at n harmonics
    
    ### Parameters:
        m_c (`float/array`): Chirp mass of each binary.
        f_orb (`float/array`): Orbital frequency of each binary at each timestep.
        ecc (`float/array`): Eccentricity of each binary at each timestep. 
        n (`int/array`): Harmonic(s) at which to calculate the strain. 
        dist (`float/array`): Distance to each binary. 
        position (`SkyCoord/array`, optional): Sky position of source. Defaults to None.
        polarisation (`float/array`, optional): GW polarisation angle of the source. Defaults to None.
        inclination (`float/array`, optional): Inclination of the source. Defaults to None.
        interpolated_g (`function`, optional): Computes g(n,e) from Peters (1964). Defaults to None.
    
    ### Returns:
        h_c (`float/array`): Strain amplitude.
    """ 
     
    arrayed_args, _ = utils.ensure_array(m_c, f_orb, ecc, n, dist) # convert to array
    m_c, f_orb, ecc, n, dist = arrayed_args

    if f_orb.ndim != 2: # extend timestep dimensions
        f_orb = f_orb[:, np.newaxis]
    if ecc.ndim != 2:
        ecc = ecc[:, np.newaxis]

    m_c = m_c[:, np.newaxis] # extend mass and distance dimensions
    dist = dist[:, np.newaxis]

    # get strain for n independent part
    prefac = (2**(5/3) / (3 * np.pi**(4/3)))**(0.5) * c.G**(5/6) / c.c**(3/2)
    n_independent_part = prefac * m_c**(5/6) / dist * f_orb**(-1/6) / utils.peters_f(ecc)**(0.5)

    if interpolated_g is None: # check whether to interpolate g(n, e)
        n = n[np.newaxis, np.newaxis, :] # extend harmonic, eccentricity dimensions
        ecc = ecc[..., np.newaxis]
        n_dependent_part = (utils.peters_g(n, ecc) / n)**(1/2)
    else:
        g_vals = interpolated_g(n, ecc.flatten()) # flatten array for interp2d
        g_vals[g_vals < 0.0] = 0.0 # set negative values from cubic fit to 0

        if isinstance(ecc, (np.ndarray, list)) and len(ecc) > 1: # unsort output array
            g_vals = g_vals[np.argsort(ecc.flatten()).argsort()] 

        g_vals = g_vals.reshape((*ecc.shape, len(n))) # reshape output dimensions
        n_dependent_part = (g_vals / n[np.newaxis, np.newaxis, :])**(0.5)

    h_c = n_independent_part[..., np.newaxis] * n_dependent_part

    if position is not None or polarisation is not None or inclination is not None:
        raise ValueError("Specific source position, polarisation, and/or inclination not available.")

    return h_c.decompose().value
