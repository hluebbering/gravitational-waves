import numpy as np
import GravitationalWaves.visualization as visualization
import unittest
from astropy import units as u


class Testvisualization(unittest.TestCase):
    """Tests that the code is functioning properly"""

    def test_approximate_response_function(self):
        """check that the R approximation works for low frequencies"""
        frequencies = np.logspace(-6, -2, 10000) * u.Hz

        exact = psd.power_spectral_density(frequencies, approximate_R = True)
        approx = psd.power_spectral_density(frequencies, approximate_R = False)

        self.assertTrue(np.allclose(exact, approx))