import numpy as np
import GravitationalWaves.psd as psd
import unittest
from astropy import units as u


class Testpsd(unittest.TestCase):
    """Tests that the code is functioning properly"""

    def test_approximate_response_function(self):
        """check that the R approximation works for low frequencies"""
        frequencies = np.logspace(-6, -2, 10000) * u.Hz

        exact = psd.power_spectral_density(frequencies, approximate_R = True)
        approx = psd.power_spectral_density(frequencies, approximate_R = False)

        self.assertTrue(np.allclose(exact, approx))

    def test_get_confusion_noise(self):
        """check that confusion noise is doing logical things"""
        frequencies = np.logspace(-6, 0, 10000) * u.Hz

        model = "robson19"
        instrument = "LISA"
        confused = psd.get_confusion_noise(frequencies, model = model)
        lucid = psd.get_confusion_noise(frequencies, model = None)

        # ensure confusion noise only adds to noise
        self.assertTrue(np.all(confused >= lucid))

        # ensure that it doesn't affect things at low or high frequency
        safe = np.logical_or(frequencies < 1e-4 * u.Hz, frequencies > 1e-2 * u.Hz)
        self.assertTrue(np.allclose(confused[safe], lucid[safe]))
        

    def test_power_spectral_density(self):
        """check that increasing the mission length isn't changing
        anything far from confusion noise"""
        frequencies = np.logspace(-6, 0, 100) * u.Hz

        # compute same curve with various mission length
        smol = psd.power_spectral_density(frequencies, t_obs=0.5 * u.yr)
        teeny = psd.power_spectral_density(frequencies, t_obs=1.0 * u.yr)
        little = psd.power_spectral_density(frequencies, t_obs=2.0 * u.yr)
        regular = psd.power_spectral_density(frequencies, t_obs=4.5 * u.yr)
        looonngg = psd.power_spectral_density(frequencies, t_obs=10.0 * u.yr)
        noises = [smol, teeny, little, regular, looonngg]

        # ensure that a shorter mission length never decreases the noise
        for noise in noises:
            above = noise > regular
            close = np.isclose(noise, regular, atol=1e-39)
            self.assertTrue(np.logical_or(above, close).all())
            
    def test_lisa_psd(self):
        """check that changing instruments doesn't break things"""
        frequencies = np.logspace(-6, 0, 100) * u.Hz

        tlp = psd.power_spectral_density(frequencies, instrument="LISA")

        def custom_instrument(f, t_obs, L, approximate_R, confusion_noise):
            return psd.lisa_psd(f, L * 2, t_obs, approximate_R, confusion_noise)

        custom = psd.power_spectral_density(frequencies, instrument="custom",
                                            custom_psd=custom_instrument,
                                            L=np.sqrt(3) * 1e5 * u.km,
                                            confusion_noise="robson19",
                                            t_obs=5 * u.yr)

        self.assertTrue(np.all(custom <= tlp))