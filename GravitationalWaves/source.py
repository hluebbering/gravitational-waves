"""Analyzes a generic set of sources that may be stationary/evolving and circular/eccentric"""  

from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
from importlib import resources
from scipy.interpolate import interp1d, interp2d
from legwork import source as scs

from GravitationalWaves import utils, psd, strain
from GravitationalWaves import visualization as vis

__all__ = ['Source']

class Source():  
    def __init__(self, m_1, m_2, ecc, dist, n_proc=1, f_orb=None, a=None, position=None, polarisation=None,
                 inclination=None, weights=None, gw_lum_tol=0.05, stat_tol=1e-2, interpolate_g=True,
                 interpolate_sc=True, sc_params={}):
        """Class for generic GW sources
        ## Parameters:
            - m_1 (`float/array`): Primary mass. Must have astropy units of mass.
            - m_2 (`float/array`): Secondary mass. Must have astropy units of mass.
            - ecc (`float/array`): Initial eccentricity
            - dist (`float/array`): Luminosity distance to source.
            - n_proc (`int`, optional): Number of processors to split eccentric evolution over. Defaults to 1.
            - f_orb (`float/array`): Orbital frequency (either a or f_orb must be supplied).
            - a (`float/array`): Semi-major axis (either a or f_orb must be supplied).
            - position (`SkyCoord/array`, optional): Sky position of source. Defaults to None.
            - polarisation (`float/array`, optional): GW polarisation angle of the source. Defaults to None.
            - inclination (`float/array`, optional): Inclination of the source. Must have astropy angular units. Defaults to None.
            - weights (`float/array`, optional): Statistical weights associated with each sample. Defaults to None.
            - gw_lum_tol (`float`, optional): Allowed error on the GW luminosity when calculating SNRs. Defaults to 0.05.
            - stat_tol (`float`, optional): Frequency change above which a binary is considered stationary. Defaults to 1e-2.
            - interpolate_g (`bool`, optional): Whether to interpolate the g(n,e) function from Peters (1964). Defaults to True.
            - interpolate_sc (`bool`, optional): Whether to interpolate the LISA sensitivity curve. Defaults to True.
            - sc_params (`dict`, optional): Parameters for interpolated sensitivity curve. 
                - Defaults: {instrument = "LISA", custom_psd = None, t_obs = "auto", L = "auto", approximate_R = False, confusion_noise = "auto"}                
        """

        # FREQUENCY (f_orb) OR SEMI-MAJOR AXIS (a) MUST BE SUPPLIED
        if f_orb is None and a is None:
            raise ValueError("Either `f_orb` or `a` must be specified")

        # RAISE ERROR OR FILL MISSING VALUES
        # If only part of the position, inclination, and polarization are supplied
        if position is None:
            if inclination is not None:
                raise ValueError(
                    "If you specify the inclination, you must also specify a sky position.")
            if polarisation is not None:
                raise ValueError(
                    "If you specify the polarisation, you must also specify a sky position.")
        else:
            if inclination is None:
                print("Generating random values for source inclinations")
                inclination = np.arcsin(
                    np.random.uniform(-1, 1, len(m_1))) * u.rad
            if polarisation is None:
                print("Generating random values for source polarisations")
                polarisation = np.random.uniform(0, 2 * np.pi, len(m_1)) * u.rad

        # POSITION MUST BE IN CORRECT COORDINATE FRAME
        if position is not None:
            if np.atleast_1d(ecc).any() > 0.0:
                raise ValueError("Only valued for circular sources")
            position = position.transform_to("heliocentrictrueecliptic")
            lon, lat, polarisation, inclination = np.atleast_1d(position.lon, position.lat, polarisation, inclination)
            position = SkyCoord(lon=lon, lat=lat, distance=dist, frame='heliocentrictrueecliptic')

        # CALCULATE ONE NOT SUPPLIED
        f_orb = utils.get_f_orb_from_a(a, m_1, m_2) if f_orb is None else f_orb
        a = utils.get_a_from_f_orb(f_orb, m_1, m_2) if a is None else a

        # ASSIGN UNITS
        # AssertionError if a parameter is missing units
        unit_args = [m_1, m_2, dist, f_orb, a]
        unit_args_str = ['m_1', 'm_2', 'dist', 'f_orb', 'a']
        for i in range(len(unit_args)):
            assert (isinstance(unit_args[i], u.quantity.Quantity)), "`{}` must have units".format(unit_args_str[i])

        # INPUTS MUST BE ARRAYS
        fixed_args, _ = utils.ensure_array(m_1, m_2, dist, f_orb, a, ecc, weights)
        m_1, m_2, dist, f_orb, a, ecc, weights = fixed_args

        # ALL ARRAY ARGUMENTS MUST BE THE SAME LENGTH
        array_args = [m_1, m_2, dist, f_orb, a, ecc]
        length_check = np.array([len(arg) != len(array_args[0]) for arg in array_args])
        if length_check.any():
            raise ValueError("All input arrays must have the same length")

        default_sc_params = {"instrument": "LISA", "custom_psd": None, "t_obs": "auto", "L": "auto", 
                             "approximate_R": False, "confusion_noise": 'auto'}
        default_sc_params.update(sc_params)
        self._sc_params = default_sc_params # Parameters for interpolated sensitivity curve. (`dict`)
        self.m_1 = m_1  
        self.m_2 = m_2  
        self.ecc = ecc  
        self.dist = dist 
        self.stat_tol = stat_tol
        self.f_orb = f_orb  
        self.a = a  
        self.n_proc = n_proc 
        self.interpolate_sc = interpolate_sc 
        self.inclination = inclination 
        self.polarisation = polarisation 
        self.weights = weights 
        self._gw_lum_tol = gw_lum_tol 
        scs.Source.set_g(interpolate_g) 
        
        self.m_c = utils.chirp_mass(m_1, m_2)  # Chirp mass. (`float/array`)
        self.n_sources = len(m_1)  # Number of sources in class (`int`)
        self.snr = None  # Signal-to-noise ratio. (`float/array`)
        self.max_snr_harmonic = None # Harmonic with the maximum sn. (`int/array`)
        self.t_merge = None
        self.merged = np.repeat(False, self.n_sources)
        scs.Source.set_sc() # Set Source sensitivity curve function: interpolation of LISA sensitivity curve using sc_params
    

    def get_snr(self, t_obs=None, instrument=None, custom_psd=None, n_step=100,
                verbose=False, re_interpolate_sc=True, which_sources=None):
        """Computes signal-to-noise ratio (SNR) for a generic binary.
        ## Parameters:
            - t_obs (`array`, optional): Observation duration. Defaults to None.
            - instrument (`string`, optional): Instrument to observe with. Default to 'LISA'.
            - custom_psd (`function`, optional): Function for computing the PSD. Defaults to None.
            - n_step (`int`, optional): Number of time steps during observation duration. Defaults to 100.
            - verbose (`bool`, optional): Whether to print additional information. Defaults to False.
            - re_interpolate_sc (`bool`, optional): Whether to re-interpolate sensitivity curve. Defaults to True.
            - which_sources (`bool/array`, optional): Mask on which sources to consider stationary. Defaults to None.
        
        ## Returns:
            - SNR (`array`): Signal-to-noise ratio
        """
        # if no values are provided, use those in sc_params
        t_obs = self._sc_params["t_obs"] if t_obs is None else t_obs
        
        if verbose:
            n_snr = len(
                which_sources[which_sources]) if which_sources is not None else self.n_sources
            print("Calculating SNR for {} sources".format(n_snr))
            print("\t{}".format(len(self.merged[self.merged])),
                  "sources have already merged")
        snr = np.zeros(self.n_sources)

        # by default calculate SNR for every source
        if which_sources is None:
            which_sources = np.repeat(True, self.n_sources)

        stat_mask = np.logical_and.reduce((scs.Source.get_source_mask(circular=None,
                                                                      stationary=True,
                                                                      t_obs=t_obs),
                                           np.logical_not(self.merged),
                                           which_sources))
        evol_mask = np.logical_and.reduce((scs.Source.get_source_mask(circular=None,
                                                                      stationary=False,
                                                                      t_obs=t_obs),
                                           np.logical_not(self.merged),
                                           which_sources))
        if stat_mask.any():
            if verbose:
                n_stat = len(snr[stat_mask])
                print("\t{} sources are stationary".format(n_stat))
            snr[stat_mask] = scs.Source.get_snr_stationary(t_obs=t_obs,
                                                           which_sources=stat_mask,
                                                           verbose=verbose)
        if evol_mask.any():
            if verbose:
                n_evol = len(snr[evol_mask])
                print("\t{} sources are evolving".format(n_evol))
            snr[evol_mask] = scs.Source.get_snr_evolving(t_obs=t_obs,
                                                         which_sources=evol_mask,
                                                         n_step=n_step,
                                                         verbose=verbose)
        return snr


    def get_h_c_n(self, harmonics, which_sources=None):
        """Computes the characteristic strain for binaries for the given harmonics.
        ## Parameters:
            - harmonics (`int/array`): Harmonic(s) at which to calculate the strain.
            - which_sources (`boolean/array`, optional): Mask on sources to compute values. Defaults to None.
        
        ## Returns:
            - h_c_n (`float/array`): Dimensionless characteristic strain.
        """        
        if which_sources is None:
            which_sources = np.repeat(True, self.n_sources)

        # by default set all strains to zero
        n_harmonics = len(harmonics) if not isinstance(harmonics, int) else 1
        h_c_n = np.zeros((self.n_sources, n_harmonics))

        # find a mask for the inspiralling sources (exclude merged)
        insp_sources = np.logical_and(
            np.logical_not(self.merged), which_sources)

        # only apply the mask to position, polarization and inclination (optional)
        position = self.position[insp_sources] if self.position is not None else None
        polarisation = self.polarisation[insp_sources] if self.position is not None else None
        inclination = self.inclination[insp_sources] if self.position is not None else None

        # calculate strain for these values
        h_c_n[insp_sources, :] = strain.h_c_n(m_c=self.m_c[insp_sources],
                                              f_orb=self.f_orb[insp_sources],
                                              ecc=self.ecc[insp_sources],
                                              n=harmonics,
                                              dist=self.dist[insp_sources],
                                              position=position,
                                              polarisation=polarisation,
                                              inclination=inclination,
                                              interpolated_g=self.g)[:, 0, :]
        # return all sources, not just inpsiralling ones
        return h_c_n[which_sources, :]

    def plot_source_variables(self, xstr, ystr=None, which_sources=None, exclude_merged_sources=True, **kwargs):
        """ Plot distributions of Source variables.
        ## Parameters:
            - xstr ({'m_1', 'm_2', 'm_c', 'ecc', 'dist', 'f_orb', 'f_GW', 'a', 'snr'}): Which variable to plot on the x axis.
            - ystr ({'m_1', 'm_2', 'm_c', 'ecc', 'dist', 'f_orb', 'f_GW', 'a', 'snr'}): Which variable to plot on the y axis.
            - which_sources (`boolean/array`, optional): Mask for which sources should be plotted. Default is all sources.
            - exclude_merged_sources (`bool`, optional): Whether to exclude merged sources in distributions. Defaults to True.
        
        ## Returns:
            - fig (`matplotlib Figure`): The figure on which the distribution is plotted.
            - ax (`matplotlib Axis`): The axis on which the distribution is plotted.
        """        
                
        convert = {"m_1": self.m_1, "m_2": self.m_2, "m_c": self.m_c,
                   "ecc": self.ecc * u.dimensionless_unscaled, "dist": self.dist, "f_orb": self.f_orb,
                   "f_GW": self.f_orb * 2, "a": self.a,
                   "snr": self.snr * u.dimensionless_unscaled if self.snr is not None else self.snr}
        labels = {"m_1": "Primary Mass", "m_2": "Secondary Mass", "m_c": "Chirp Mass", "ecc": "Eccentricity",
                  "dist": "Distance", "f_orb": "Orbital Frequency", "f_GW": "Gravitational Wave Frequency",
                  "a": "Semi-major axis", "snr": "Signal-to-noise Ratio"}
        unitless = set(["ecc", "snr"])

        if which_sources is None:
            which_sources = np.repeat(True, self.n_sources)

        if exclude_merged_sources:
            which_sources = np.logical_and(which_sources, np.logical_not(self.merged))

        # ensure that the variable is a valid choice
        for var_str in [xstr, ystr]:
            if var_str not in convert.keys() and var_str is not None:
                error_str = "`xstr` and `ystr` must be one of: " + ', '.join(["`{}`".format(k) for k in list(convert.keys())])
                raise ValueError(error_str)

        # check the instance variable has been already set
        x = convert[xstr]
        if x is None:
            raise ValueError("x variable (`{}`)".format(xstr), "must be not be None")
        if ystr is not None:
            y = convert[ystr]
            if y is None:
                raise ValueError("y variable (`{}`)".format(ystr), "must be not be None")

        # create the x label if it wasn't provided
        if "xlabel" not in kwargs.keys():
            if xstr in unitless:
                kwargs["xlabel"] = labels[xstr]
            else:
                kwargs["xlabel"] = r"{} [{:latex}]".format(labels[xstr], x.unit)

        # create the y label if it wasn't provided and ystr was
        if ystr is not None and "ylabel" not in kwargs.keys():
            if ystr in unitless:
                kwargs["ylabel"] = labels[ystr]
            else:
                kwargs["ylabel"] = r"{} [{:latex}]".format(labels[ystr], y.unit)

        # work out what the weights are
        weights = self.weights[which_sources] if self.weights is not None else None

        # If two variables are specified then produce a 2D distribution
        if ystr is not None:
            return vis.plot_2D_dist(x=x[which_sources], y=y[which_sources], weights=weights, **kwargs)
        else:
            return vis.plot_1D_dist(x=x[which_sources], weights=weights, **kwargs)
