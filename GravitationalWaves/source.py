from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
from importlib import resources
from scipy.interpolate import interp1d, interp2d
from legwork import evol

from GravitationalWaves import utils, psd, strain
import GravitationalWaves.snr as sn
import GravitationalWaves.visualization as vis

__all__ = ['Source']


# Class for generic GW sources
# Analyzes a generic set of sources that may be stationary/evolving and circular/eccentric
class Source():
    def __init__(self, m_1, m_2, ecc, dist, n_proc=1, f_orb=None, a=None, position=None, polarisation=None,
                 inclination=None, weights=None, gw_lum_tol=0.05, stat_tol=1e-2, interpolate_g=True,
                 interpolate_sc=True, sc_params={}):

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
                polarisation = np.random.uniform(
                    0, 2 * np.pi, len(m_1)) * u.rad

        # POSITION MUST BE IN CORRECT COORDINATE FRAME
        if position is not None:
            if np.atleast_1d(ecc).any() > 0.0:
                raise ValueError("Only valued for circular sources")
            position = position.transform_to("heliocentrictrueecliptic")
            # POSITION, POLARISATION, AND INCLINATION ARE 1D
            lon, lat, polarisation, inclination = np.atleast_1d(
                position.lon, position.lat, polarisation, inclination)
            position = SkyCoord(lon=lon, lat=lat, distance=dist,
                                frame='heliocentrictrueecliptic')

        # CALCULATE ONE NOT SUPPLIED
        f_orb = utils.get_f_orb_from_a(a, m_1, m_2) if f_orb is None else f_orb
        a = utils.get_a_from_f_orb(f_orb, m_1, m_2) if a is None else a

        # ASSIGN UNITS
        # AssertionError if a parameter is missing units
        unit_args = [m_1, m_2, dist, f_orb, a]
        unit_args_str = ['m_1', 'm_2', 'dist', 'f_orb', 'a']
        for i in range(len(unit_args)):
            assert (isinstance(unit_args[i], u.quantity.Quantity)), \
                "`{}` must have units".format(unit_args_str[i])

        # INPUTS MUST BE ARRAYS
        fixed_args, _ = utils.ensure_array(
            m_1, m_2, dist, f_orb, a, ecc, weights)
        m_1, m_2, dist, f_orb, a, ecc, weights = fixed_args

        # ALL ARRAY ARGUMENTS MUST BE THE SAME LENGTH
        array_args = [m_1, m_2, dist, f_orb, a, ecc]
        length_check = np.array([len(arg) != len(array_args[0])
                                for arg in array_args])
        if length_check.any():
            raise ValueError("All input arrays must have the same length")

        default_sc_params = {
            "instrument": "LISA",
            "custom_psd": None,
            "t_obs": "auto",
            "L": "auto",
            "approximate_R": False,
            "confusion_noise": 'auto'
        }
        default_sc_params.update(sc_params)

        # PARAMETERS
        self._sc_params = default_sc_params # Parameters for interpolated sensitivity curve. (`dict`)
        self.m_1 = m_1  # Primary mass. (`float/array`)
        self.m_2 = m_2  # Secondary mass. (`float/array`)
        self.ecc = ecc  # Initial eccentricity (`float/array`)
        self.dist = dist  # Luminosity distance to source. (`float/array`)
        # Change in frequency above which a binary should be considered stationary (`float`)
        self.stat_tol = stat_tol
        self.f_orb = f_orb  # Orbital frequency (`float/array`)
        self.a = a  # Semi-major axis (`float/array`)
        self.n_proc = n_proc # Number of processors to split eccentric evolution over (`int`)
        self.interpolate_sc = interpolate_sc # Whether to interpolate the LISA sensitivity curve (`boolean`)
        self.inclination = inclination # Inclination of the source. (`float/array`)
        self.polarisation = polarisation # GW polarisation angle of source. (`float/array`, optional)
        self.weights = weights # Statistical weights associated with each sample (`float/array`, optional)
        self.update_gw_lum_tol(gw_lum_tol) # Allowed error on the GW luminosity when calculating SNRs. (`float`)
        self.set_g(interpolate_g) # Whether to interpolate the g(n,e) function from Peters (1964). (`boolean`)

        # ATTRIBUTES
        self.m_c = utils.chirp_mass(m_1, m_2)  # Chirp mass. (`float/array`)
        self.n_sources = len(m_1)  # Number of sources in class (`int`)
        self.snr = None  # Signal-to-noise ratio. (`float/array`)
        self.max_snr_harmonic = None # Harmonic with the maximum sn. (`int/array`)
        self.t_merge = None
        self.merged = np.repeat(False, self.n_sources)
        self.set_sc()

    # Set Source g function if user wants to interpolate g(n,e)
    # Parameters: interpolate_g (boolean) - Whether to interpolate the g(n,e) function from Peters (1964)
    def set_g(self, interpolate_g):
        if interpolate_g:
            # open file with pre-calculated g(n,e) grid
            with resources.path(package="legwork", resource="peters_g.npy") as path:
                peters_g = np.load(path)
            # interpolate grid using scipy
            n_range = np.arange(1, 10000 + 1).astype(int)
            e_range = np.linspace(0, 1, 1000)
            self.g = interp2d(n_range, e_range, peters_g, kind="cubic")
        else:
            self.g = None

    # Update GW luminosity tolerance to recalculate harmonics_required function and transition to eccentric
    # Parameters: gw_lum_tol (float) - Allowed error on the GW luminosity when calculating SNRs
    def update_gw_lum_tol(self, gw_lum_tol):
        self._gw_lum_tol = gw_lum_tol
        self.create_harmonics_functions()
        self.find_eccentric_transition()

    # Set Source sensitivity curve function: interpolation of LISA sensitivity curve using sc_params
    def set_sc(self):
        if self.interpolate_sc:
            # get values
            frequency_range = np.logspace(-7, np.log10(2), 10000) * u.Hz
            sc = psd.power_spectral_density(frequency_range, **self._sc_params)
            # interpolate
            interp_sc = interp1d(frequency_range, sc,
                                 bounds_error=False, fill_value=1e30)
            # add units back
            self.sc = lambda f: interp_sc(f.to(u.Hz)) / u.Hz
        else:
            self.sc = None

    # Produce a mask of the sources based on whether binaries are circular or eccentric and stationary or evolving
    def get_source_mask(self, circular=None, stationary=None, t_obs=None):
        """Parameters: 
        circular (bool): None = either, True = only circular binaries, False = only eccentric
        stationary (bool): None = either, True = only stationary binaries, False = only evolving
        t_obs (float): Observation time
        """

        if circular is None:
            circular_mask = np.repeat(True, self.n_sources)
        elif circular is True:
            circular_mask = self.ecc <= self.ecc_tol
        elif circular is False:
            circular_mask = self.ecc > self.ecc_tol
        else:
            raise ValueError("`circular` must be None, True or False")

        if stationary is None:
            stat_mask = np.repeat(True, self.n_sources)
        elif stationary is True or stationary is False:
            t_obs = self._sc_params["t_obs"] if t_obs is None else t_obs
            if t_obs == "auto":
                t_obs = 4 * \
                    u.yr if self._sc_params["instrument"] == "LISA" else 5 * u.yr
            stat_mask = evol.determine_stationarity(m_c=self.m_c, f_orb_i=self.f_orb, t_evol=t_obs,
                                                    ecc_i=self.ecc, stat_tol=self.stat_tol)
            if stationary is False:
                stat_mask = np.logical_not(stat_mask)
        else:
            raise ValueError("`stationary` must be None, True or False")

        # Returns: mask (bool/array)
        return np.logical_and(circular_mask, stat_mask)

    # Computes the SNR assuming a stationary binary
    def get_snr_stationary(self, t_obs=None, instrument=None, custom_psd=None, which_sources=None, verbose=False):
        """Parameters: 
        t_obs (array): Observation duration
        instrument ('LISA', 'TianQin', 'custom') 
        custom_psd (function): Function for computing the PSD
        which_sources (bool/array): Mask on which sources to consider stationary
        verbose (boolean): Whether to print additional information
        """
        if which_sources is None:
            which_sources = np.repeat(True, self.n_sources)

        instrument = self._sc_params["instrument"] if instrument is None else instrument
        t_obs = self._sc_params["t_obs"] if t_obs is None else t_obs
        if t_obs == "auto":
            t_obs = 4 * u.yr if instrument == "LISA" else 5 * u.yr

        insp_sources = np.logical_and(
            which_sources, np.logical_not(self.merged))
        snr = np.zeros(self.n_sources)
        e_mask = np.logical_and(self.ecc > self.ecc_tol, insp_sources)
        c_mask = np.logical_and(self.ecc <= self.ecc_tol, insp_sources)

        # default to n = 2 for max snr harmonic
        msh = np.repeat(2, self.n_sources)

        # only apply the mask to the position, polarization and inclination
        # if they are provided
        position = self.position[c_mask] if self.position is not None else None
        polarisation = self.polarisation[c_mask] if self.position is not None else None
        inclination = self.inclination[c_mask] if self.position is not None else None

        # only compute snr if there is at least one binary in mask
        if c_mask.any():
            if verbose:
                print("\t\t{} sources are stationary and circular".format(
                    len(snr[c_mask])))
            snr[c_mask] = sn.snr_circ_stationary(m_c=self.m_c[c_mask],
                                                 f_orb=self.f_orb[c_mask],
                                                 dist=self.dist[c_mask],
                                                 t_obs=t_obs,
                                                 interpolated_g=self.g,
                                                 interpolated_sc=self.sc,
                                                 instrument=instrument,
                                                 custom_psd=custom_psd,
                                                 position=position,
                                                 polarisation=polarisation,
                                                 inclination=inclination)
        if e_mask.any():
            if verbose:
                print("\t\t{} sources are stationary and eccentric".format(
                    len(snr[e_mask])))
            harmonics_required = self.harmonics_required(self.ecc)
            harmonic_groups = [(1, 10), (10, 100), (100, 1000), (1000, 10000)]
            for lower, upper in harmonic_groups:
                harm_mask = np.logical_and(
                    harmonics_required > lower, harmonics_required <= upper)
                match = np.logical_and(harm_mask, e_mask)

                # only apply the mask to the position, polarization and inclination
                # if they are provided
                position = self.position[match] if self.position is not None else None
                polarisation = self.polarisation[match] if self.position is not None else None
                inclination = self.inclination[match] if self.position is not None else None

                if match.any():
                    snr[match], msh[match] = sn.snr_ecc_stationary(m_c=self.m_c[match],
                                                                   f_orb=self.f_orb[match],
                                                                   ecc=self.ecc[match],
                                                                   dist=self.dist[match],
                                                                   t_obs=t_obs,
                                                                   harmonics_required=upper,
                                                                   interpolated_g=self.g,
                                                                   interpolated_sc=self.sc,
                                                                   ret_max_snr_harmonic=True,
                                                                   instrument=instrument,
                                                                   custom_psd=custom_psd,
                                                                   position=position,
                                                                   polarisation=polarisation,
                                                                   inclination=inclination)

        if self.max_snr_harmonic is None:
            self.max_snr_harmonic = np.zeros(self.n_sources).astype(int)
        self.max_snr_harmonic[insp_sources] = msh[insp_sources]

        if self.snr is None:
            self.snr = np.zeros(self.n_sources)
        self.snr[insp_sources] = snr[insp_sources]

        # Returns: SNR (array): Signal-to-noise ratio
        return snr[which_sources]

    # Computes the SNR assuming an evolving binary
    def get_snr_evolving(self, t_obs=None, instrument=None, custom_psd=None, n_step=100, which_sources=None,
                         verbose=False):
        """Parameters: 
        t_obs (array): Observation duration
        instrument ('LISA', 'TianQin', 'custom') 
        custom_psd (function): Function for computing the PSD
        which_sources (bool/array): Mask on which sources to consider stationary
        verbose (boolean): Whether to print additional information
        """
        snr = np.zeros(self.n_sources)

        if which_sources is None:
            which_sources = np.repeat(True, self.n_sources)

        instrument = self._sc_params["instrument"] if instrument is None else instrument
        t_obs = self._sc_params["t_obs"] if t_obs is None else t_obs
        if t_obs == "auto":
            t_obs = 4 * u.yr if instrument == "LISA" else 5 * u.yr

        insp_sources = np.logical_and(
            which_sources, np.logical_not(self.merged))
        e_mask = np.logical_and(self.ecc > self.ecc_tol, insp_sources)
        c_mask = np.logical_and(self.ecc <= self.ecc_tol, insp_sources)

        # only apply the mask to the position, polarization and inclination
        # if they are provided
        position = self.position[c_mask] if self.position is not None else None
        polarisation = self.polarisation[c_mask] if self.position is not None else None
        inclination = self.inclination[c_mask] if self.position is not None else None

        # default to n = 2 for max snr harmonic
        msh = np.repeat(2, self.n_sources)

        if c_mask.any():
            if verbose:
                print("\t\t{} sources are evolving and circular".format(
                    len(snr[c_mask])))
            t_merge = None if self.t_merge is None else self.t_merge[c_mask]
            snr[c_mask] = sn.snr_circ_evolving(m_1=self.m_1[c_mask],
                                               m_2=self.m_2[c_mask],
                                               f_orb_i=self.f_orb[c_mask],
                                               dist=self.dist[c_mask],
                                               t_obs=t_obs,
                                               n_step=n_step,
                                               t_merge=t_merge,
                                               interpolated_g=self.g,
                                               interpolated_sc=self.sc,
                                               instrument=instrument,
                                               custom_psd=custom_psd,
                                               position=position,
                                               polarisation=polarisation,
                                               inclination=inclination)
        if e_mask.any():
            if verbose:
                print("\t\t{} sources are evolving and eccentric".format(
                    len(snr[e_mask])))
            harmonics_required = self.harmonics_required(self.ecc)
            harmonic_groups = [(1, 10), (10, 100), (100, 1000), (1000, 10000)]
            for lower, upper in harmonic_groups:
                harm_mask = np.logical_and(
                    harmonics_required > lower, harmonics_required <= upper)
                match = np.logical_and(harm_mask, e_mask)
                # only apply the mask to the position, polarization and inclination
                # if they are provided
                position = self.position[match] if self.position is not None else None
                polarisation = self.polarisation[match] if self.position is not None else None
                inclination = self.inclination[match] if self.position is not None else None

                if match.any():
                    t_merge = None if self.t_merge is None else self.t_merge[match]
                    snr[match], msh[match] = sn.snr_ecc_evolving(m_1=self.m_1[match],
                                                                 m_2=self.m_2[match],
                                                                 f_orb_i=self.f_orb[match],
                                                                 dist=self.dist[match],
                                                                 ecc=self.ecc[match],
                                                                 harmonics_required=upper,
                                                                 t_obs=t_obs,
                                                                 n_step=n_step,
                                                                 interpolated_g=self.g,
                                                                 interpolated_sc=self.sc,
                                                                 n_proc=self.n_proc,
                                                                 ret_max_snr_harmonic=True,
                                                                 instrument=instrument,
                                                                 custom_psd=custom_psd,
                                                                 position=position,
                                                                 polarisation=polarisation,
                                                                 inclination=inclination)

        if self.max_snr_harmonic is None:
            self.max_snr_harmonic = np.zeros(self.n_sources).astype(int)
        self.max_snr_harmonic[insp_sources] = msh[insp_sources]

        if self.snr is None:
            self.snr = np.zeros(self.n_sources)
        self.snr[insp_sources] = snr[insp_sources]

        # Returns: SNR (array): Signal-to-noise ratio
        return snr[which_sources]

    # Computes SNR for a generic binary.
    def get_snr(self, t_obs=None, instrument=None, custom_psd=None, n_step=100,
                verbose=False, re_interpolate_sc=True, which_sources=None):
        """Parameters: 
        t_obs (array): Observation duration
        instrument ('LISA', 'TianQin', 'custom') 
        custom_psd (function): Function for computing the PSD
        which_sources (bool/array): Mask on which sources to consider stationary
        verbose (boolean): Whether to print additional information
        """
        # if no values are provided, use those in sc_params
        t_obs = self._sc_params["t_obs"] if t_obs is None else t_obs
        instrument = self._sc_params["instrument"] if instrument is None else instrument
        custom_psd = self._sc_params["custom_psd"] if custom_psd is None else custom_psd

        # if the user interpolated a sensitivity curve with different settings
        if (self.interpolate_sc and self._sc_params is not None
                and (t_obs != self._sc_params["t_obs"]
                     or instrument != self._sc_params["instrument"]
                     or custom_psd != self._sc_params["custom_psd"])):  # pragma: no cover

            # re interpolate the sensitivity curve with new parameters
            if re_interpolate_sc:
                self._sc_params["t_obs"] = t_obs
                self._sc_params["instrument"] = instrument
                self._sc_params["custom_psd"] = custom_psd

                self.set_sc()

            # otherwise warn the user that they are making a mistake
            else:
                print("Set `re_interpolate_sc=True` to re-interpolate the sensitivity curve.")

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

        stat_mask = np.logical_and.reduce((self.get_source_mask(circular=None,
                                                                stationary=True,
                                                                t_obs=t_obs),
                                           np.logical_not(self.merged),
                                           which_sources))
        evol_mask = np.logical_and.reduce((self.get_source_mask(circular=None,
                                                                stationary=False,
                                                                t_obs=t_obs),
                                           np.logical_not(self.merged),
                                           which_sources))

        if stat_mask.any():
            if verbose:
                n_stat = len(snr[stat_mask])
                print("\t{} sources are stationary".format(n_stat))
            snr[stat_mask] = self.get_snr_stationary(t_obs=t_obs,
                                                     instrument=instrument,
                                                     custom_psd=custom_psd,
                                                     which_sources=stat_mask,
                                                     verbose=verbose)
        if evol_mask.any():
            if verbose:
                n_evol = len(snr[evol_mask])
                print("\t{} sources are evolving".format(n_evol))
            snr[evol_mask] = self.get_snr_evolving(t_obs=t_obs,
                                                   instrument=instrument,
                                                   custom_psd=custom_psd,
                                                   which_sources=evol_mask,
                                                   n_step=n_step,
                                                   verbose=verbose)
        # Returns: SNR (array): Signal-to-noise ratio
        return snr

    # Computes the characteristic strain for binaries for the given harmonics.
    def get_h_c_n(self, harmonics, which_sources=None):
        """Parameters:
        harmonics (int/array): Harmonic(s) at which to calculate the strain
        which_sources (boolean/array): Mask on sources to compute values
        """
        if which_sources is None:
            which_sources = np.repeat(True, self.n_sources)

        # by default set all strains to zero
        n_harmonics = len(harmonics) if not isinstance(harmonics, int) else 1
        h_c_n = np.zeros((self.n_sources, n_harmonics))

        # find a mask for the inspiralling sources (exclude merged)
        insp_sources = np.logical_and(
            np.logical_not(self.merged), which_sources)

        # only apply the mask to the position, polarization and inclination
        # if they are provided
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

        # Returns: h_c_n (float/array): Dimensionless characteristic strain
        # return all sources, not just inpsiralling ones
        return h_c_n[which_sources, :]

    # Plot distributions of Source variables.
    def plot_source_variables(self, xstr, ystr=None, which_sources=None,
                              exclude_merged_sources=True, **kwargs):  # pragma: no cover
        """Parameters
        xstr : `{ 'm_1', 'm_2', 'm_c', 'ecc', 'dist', 'f_orb', 'f_GW', 'a', 'snr' }`
            Which variable to plot on the x axis
        ystr : `{ 'm_1', 'm_2', 'm_c', 'ecc', 'dist', 'f_orb', 'f_GW', 'a', snr' }`
            Which variable to plot on the y axis (if None then a 1D distribution is made using `xstr`)
        which_sources : `boolean array`
            Mask for which sources should be plotted
        exclude_merged_sources : `boolean`
            Whether to exclude merged sources in distributions
        **kwargs : `various`
            When only ``xstr`` is provided, the kwargs are the same as
            :meth:`legwork.visualisation.plot_1D_dist`. When both ``xstr`` and ``ystr`` are provided,
            the kwargs are the same as :meth:`legwork.visualisation.plot_2D_dist`.
            Note that if ``xlabel`` or ``ylabel`` is not passed then this function automatically creates
            one using a default string and (if applicable) the Astropy units of the variable.
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
            which_sources = np.logical_and(
                which_sources, np.logical_not(self.merged))

        # ensure that the variable is a valid choice
        for var_str in [xstr, ystr]:
            if var_str not in convert.keys() and var_str is not None:
                error_str = "`xstr` and `ystr` must be one of: " \
                            + ', '.join(["`{}`".format(k)
                                        for k in list(convert.keys())])
                raise ValueError(error_str)

        # check the instance variable has been already set
        x = convert[xstr]
        if x is None:
            raise ValueError("x variable (`{}`)".format(
                xstr), "must be not be None")
        if ystr is not None:
            y = convert[ystr]
            if y is None:
                raise ValueError("y variable (`{}`)".format(
                    ystr), "must be not be None")

        # create the x label if it wasn't provided
        if "xlabel" not in kwargs.keys():
            if xstr in unitless:
                kwargs["xlabel"] = labels[xstr]
            else:
                kwargs["xlabel"] = r"{} [{:latex}]".format(
                    labels[xstr], x.unit)

        # create the y label if it wasn't provided and ystr was
        if ystr is not None and "ylabel" not in kwargs.keys():
            if ystr in unitless:
                kwargs["ylabel"] = labels[ystr]
            else:
                kwargs["ylabel"] = r"{} [{:latex}]".format(
                    labels[ystr], y.unit)

        # work out what the weights are
        weights = self.weights[which_sources] if self.weights is not None else None

        # Returns:
        # fig (matplotlib Figure)
        #   The figure on which the distribution is plotted
        # ax (matplotlib Axis)
        #   The axis on which the distribution is plotted

        # If two variables are specified then produce a 2D distribution
        if ystr is not None:
            return vis.plot_2D_dist(x=x[which_sources], y=y[which_sources], weights=weights, **kwargs)
        # Otherwise a 1D distribution
        else:
            return vis.plot_1D_dist(x=x[which_sources], weights=weights, **kwargs)
