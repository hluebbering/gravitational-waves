import GravitationalWaves.source as source
import GravitationalWaves.strain as strain
import GravitationalWaves.utils as utils
import numpy as np
import legwork.snr as snr
import unittest

from astropy import units as u
from astropy.coordinates import SkyCoord


class Testsource(unittest.TestCase):
    """Tests that the code is functioning properly"""

    def test_source_snr(self):
        """check that source calculates snr in correct way"""

        # create random (circular/stationary) binaries
        n_values = 500
        t_obs = 4 * u.yr
        m_1 = np.random.uniform(0, 10, n_values) * u.Msun
        m_2 = np.random.uniform(0, 10, n_values) * u.Msun
        m_c = utils.chirp_mass(m_1, m_2)
        dist = np.random.uniform(0, 30, n_values) * u.kpc
        f_orb = 10**(np.random.uniform(-5, -4, n_values)) * u.Hz
        ecc = np.repeat(0.0, n_values)
        sources = source.Source(m_1=m_1, m_2=m_2, f_orb=f_orb,
                                ecc=ecc, dist=dist)

        # compare snr calculated directly with through Source
        snr_direct = snr.snr_circ_stationary(m_c=m_c, f_orb=f_orb,
                                             dist=dist, t_obs=t_obs)
        snr_source = sources.get_snr(verbose=True)

        self.assertTrue(np.allclose(snr_direct, snr_source))

        # repeat the same test for eccentric systems
        ecc = np.random.uniform(sources.ecc_tol, 0.1, n_values)
        sources.ecc = ecc

        snr_direct = snr.snr_ecc_stationary(m_c=m_c, f_orb=f_orb, ecc=ecc,
                                            dist=dist, t_obs=t_obs,
                                            harmonics_required=10)
        snr_source = sources.get_snr(verbose=True)

        self.assertTrue(np.allclose(snr_direct, snr_source))

    def test_source_strain(self):
        """check that source calculate strain correctly"""
        n_values = 500
        m_1 = np.random.uniform(0, 10, n_values) * u.Msun
        m_2 = np.random.uniform(0, 10, n_values) * u.Msun
        m_c = utils.chirp_mass(m_1, m_2)
        dist = np.random.uniform(0, 30, n_values) * u.kpc
        f_orb = 10**(np.random.uniform(-5, -4, n_values)) * u.Hz
        ecc = np.repeat(0.0, n_values)

        sources = source.Source(m_1=m_1, m_2=m_2, f_orb=f_orb,
                                ecc=ecc, dist=dist, interpolate_g=False)

        source_strain = sources.get_h_0_n([1, 2, 3])
        true_strain = strain.h_0_n(m_c=m_c, f_orb=f_orb, ecc=ecc,
                                   n=[1, 2, 3], dist=dist)[:, 0, :]

        self.assertTrue(np.all(source_strain == true_strain))

        source_char_strain = sources.get_h_c_n([1, 2, 3])
        true_char_strain = strain.h_c_n(m_c=m_c, f_orb=f_orb, ecc=ecc,
                                        n=[1, 2, 3], dist=dist)[:, 0, :]

        self.assertTrue(np.all(source_char_strain == true_char_strain))
        
        
    def test_amplitude_modulation_h_c_n(self):
        """Make sure that the amplitude modulated strains are correct.
        Note that this is very redundant with the utils modulation tests"""
        n_values = 500
        m_1 = np.random.uniform(0, 10, n_values) * u.Msun
        m_2 = np.random.uniform(0, 10, n_values) * u.Msun
        m_c = utils.chirp_mass(m_1, m_2)
        dist = np.random.uniform(0, 30, n_values) * u.kpc
        f_orb = 10**(np.random.uniform(-5, -4, n_values)) * u.Hz
        ecc = np.repeat(0.0, n_values)
        incs = np.arccos(np.random.uniform(-1, 1, n_values)) * u.rad
        thetas = np.arcsin(np.random.uniform(-1, 1, n_values)) * u.rad
        phis = np.random.uniform(0, 2 * np.pi, n_values) * u.rad
        psis = np.random.uniform(0, 2 * np.pi, n_values) * u.rad

        positions = SkyCoord(phis, thetas, distance=dist, frame='heliocentrictrueecliptic')

        sources = source.Source(m_1=m_1, m_2=m_2, f_orb=f_orb,
                                ecc=ecc, dist=dist,
                                position=positions, inclination=incs,
                                polarisation=psis, interpolate_g=False)
        source_strains = sources.get_h_c_n([1, 2, 3])
        true_strain = strain.h_c_n(m_c=m_c, f_orb=f_orb, ecc=ecc, dist=dist, position=positions,
                                   inclination=incs, polarisation=psis, n=[1, 2, 3])[:, 0, :]
        self.assertTrue(np.all(source_strains == true_strain))
    