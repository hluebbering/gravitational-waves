import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import astropy.units as u
import legwork.psd as psd
from astropy.visualization import quantity_support

# set the default font and fontsize
plt.rc('font', family='serif')
plt.rcParams['text.usetex'] = False
fs = 24

# update various fontsizes to match
params = {'figure.figsize': (12, 8),
          'legend.fontsize': fs,
          'axes.labelsize': fs,
          'xtick.labelsize': 0.7 * fs,
          'ytick.labelsize': 0.7 * fs}
plt.rcParams.update(params)

__all__ = ['plot_1D_dist', 'plot_2D_dist', 'plot_sensitivity_curve',
           'plot_sources_on_sc_circ_stat', 'plot_sources_on_sc_ecc_stat']


def plot_1D_dist(x, weights=None, disttype="hist", fig=None, ax=None, xlabel=None, ylabel=None,
                 xlim=None, ylim=None, color=None, show=True, **kwargs):
    """plot a 1D distribution of ``x``.
    This function is a wrapper for :func:`matplotlib.pyplot.hist`, :func:`seaborn.kdeplot`
    and :func:`seaborn.ecdfplot`.
    Parameters
    ----------
    x : `float/int array`
        Variable to plot, should be a 1D array
    weights : `float/int array`
        Weights for each variable in ``x``, must have the same shape
    disttype : `{{ "hist", "kde", "ecdf" }}`
        Which type of distribution plot to use
    fig: `matplotlib Figure`
        A figure on which to plot the distribution. Both `ax` and `fig` must be supplied for either to be used
    ax: `matplotlib Axis`
        An axis on which to plot the distribution. Both `ax` and `fig` must be supplied for either to be used
    xlabel : `string`
        Label for the x axis, passed to Axes.set_xlabel()
    ylabel : `string`
        Label for the y axis, passed to Axes.set_ylabel()
    xlim : `tuple`
        Lower and upper limits for the x axis, passed to Axes.set_xlim()
    ylim : `tuple`
        Lower and upper limits for the y axis, passed to Axes.set_ylim()
    color : `string or tuple`
        Colour to use for the plot, see https://matplotlib.org/tutorials/colors/colors.html for details on
        how to specify a colour
    show : `boolean`
        Whether to immediately show the plot or only return the Figure and Axis
    **kwargs : `(if disttype=="hist")`
        Include values for any of `bins, range, density, cumulative, bottom, histtype, align, orientation,
        rwidth, log, label` or more. See :func:`matplotlib.pyplot.hist` for more details.
    **kwargs : `(if disttype=="kde")`
        Include values for any of `gridsize, cut, clip, legend, cumulative, bw_method, bw_adjust, log_scale,
        fill, label, linewidth, linestyle` or more. See :func:`seaborn.kdeplot` for more details.
    **kwargs : `(if disttype=="ecdf")`
        Include values for any of `stat, complementary, log_scale, legend, label, linewidth, linestyle` or
        more. See :func:`seaborn.edcfplot` for more details.
    Returns
    -------
    fig : `matplotlib Figure`
        The figure on which the distribution is plotted
    ax : `matplotlib Axis`
        The axis on which the distribution is plotted
    """
    # create new figure and axes is either weren't provided
    if fig is None or ax is None:
        fig, ax = plt.subplots()

    # change default kwargs for matplotlib.hist
    hist_args = {"bins": "auto", "density": True}

    # change default kwargs for seaborn.kdeplot
    kde_args = {"gridsize": 200, "legend": True}

    # change default kwargs for seaborn.ecdfplot
    ecdf_args = {"stat": 'proportion', "legend": True}

    # set which defaults we are using for this plot
    plot_args = hist_args if disttype == "hist" else kde_args if disttype == "kde" else ecdf_args

    # update the defaults with those supplied
    for key, value in kwargs.items():
        plot_args[key] = value

    # create whichever plot was requested
    with quantity_support():
        if disttype == "hist":
                ax.hist(x, weights=weights, color=color, **plot_args)
        elif disttype == "kde":
            sns.kdeplot(x=x, weights=weights, ax=ax, color=color, **plot_args)
        elif disttype == "ecdf":
            sns.ecdfplot(x=x, weights=weights, ax=ax, color=color, **plot_args)

    # update axis labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    # update axis limits
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    # immediately show the plot if requested
    if show:
        plt.show()

    # return the figure and axis for further plotting
    return fig, ax



def plot_2D_dist(x, y, weights=None, disttype="scatter", scatter_s=20, fig=None, ax=None,
                 xlabel=None, ylabel=None, xlim=None, ylim=None, color=None,
                 show=True, **kwargs):
    """Plot a 2D distribution of `x` and `y`
    This function is a wrapper for :func:`matplotlib.pyplot.scatter` and :func:`seaborn.kdeplot`.
    Parameters
    ----------
    x : `float/int array`
        Variable to plot on the x axis, should be a 1D array
    y : `float/int array`
        Variable to plot on the y axis, should be a 1D array
    weights : `float/int array`
        Weights for each variable pair (``x``, ``y``), must have the same shape
    disttype : `{{ "scatter", "kde" }}`
        Which type of distribution plot to use
    scatter_s : `float`, default=20
        Scatter point size, passed as ``s`` to a scatter plot and ignored for a KDE
    fig: `matplotlib Figure`
        A figure on which to plot the distribution. Both `ax` and `fig` must be supplied for either to be used
    ax: `matplotlib Axis`
        An axis on which to plot the distribution. Both `ax` and `fig` must be supplied for either to be used
    xlabel : `string`
        Label for the x axis, passed to Axes.set_xlabel()
    ylabel : `string`
        Label for the y axis, passed to Axes.set_ylabel()
    xlim : `tuple`
        Lower and upper limits for the x axis, passed to Axes.set_xlim()
    ylim : `tuple`
        Lower and upper limits for the u axis, passed to Axes.set_ylim()
    color : `string or tuple`
        Colour to use for the plot, see https://matplotlib.org/tutorials/colors/colors.html for details on how
        to specify a colour
    show : `boolean`
        Whether to immediately show the plot or only return the Figure and Axis
    **kwargs : `(if disttype=="scatter")`
        Input any of `s, c, marker, cmap, norm, vmin, vmax, alpha, linewidths, edgecolors` or more.
        See :func:`matplotlib.pyplot.scatter` for more details.
    **kwargs : `(if disttype=="kde")`
        Input any of `gridsize, cut, clip, legend, cumulative, cbar, cbar_ax, cbar_kws, bw_method, hue,
        palette, hue_order, hue_norm, levels, thresh, bw_adjust, log_scale, fill, label`.
        See :func:`seaborn.kdeplot` for more details.
    Returns
    -------
    fig : `matplotlib Figure`
        The figure on which the distribution is plotted
    ax : `matplotlib Axis`
        The axis on which the distribution is plotted
    """
    # create new figure and axes is either weren't provided
    if fig is None or ax is None:
        fig, ax = plt.subplots()

    # change default kwargs for matplotlib.scatter
    scatter_args = {}

    # possible kwargs for seaborn.kdeplot
    kde_args = {"gridsize": 200, "legend": True, "levels": 10, "thresh": 0.02}

    # set which ones we are using for this plot
    plot_args = scatter_args if disttype == "scatter" else kde_args

    # update the values with those supplied
    for key, value in kwargs.items():
        plot_args[key] = value

    # create whichever plot was requested
    with quantity_support():
        if disttype == "scatter":
            # change the size of scatter points based on their weights
            if weights is not None:
                scatter_s = weights * scatter_s

            ax.scatter(x, y, s=scatter_s, color=color, **plot_args)
        elif disttype == "kde":
            sns.kdeplot(x=x, y=y, weights=weights, ax=ax, color=color, **plot_args)

    # update axis labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    # update axis limits
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    # immediately show the plot if requested
    if show:
        plt.show()

    # return the figure and axis for further plotting
    return fig, ax


def plot_sensitivity_curve(frequency_range=None, y_quantity="ASD", fig=None, ax=None, show=True,
                           color="#18068b", fill=True, alpha=0.2, linewidth=1, label=None, **kwargs):
    """Plot the LISA sensitivity curve
    Parameters
    ----------
    frequency_range : `float array`
        Frequency values at which to plot the sensitivity curve
    y_quantity : `{{ "ASD", "h_c" }}`
        Which quantity to plot on the y axis (amplitude spectral density or characteristic strain)
    fig: `matplotlib Figure`
        A figure on which to plot the distribution. Both `ax` and `fig` must be supplied for either to be used
    ax: `matplotlib Axis`
        An axis on which to plot the distribution. Both `ax` and `fig` must be supplied for either to be used
    show : `boolean`
        Whether to immediately show the plot or only return the Figure and Axis
    color : `string or tuple`
        Colour to use for the curve, see https://matplotlib.org/tutorials/colors/colors.html for details on
        how to specify a colour
    fill : `boolean`
        Whether to fill the area below the sensitivity curve
    alpha : `float`
        Opacity of the filled area below the sensitivity curve (ignored if fill is `False`)
    linewidth : `float`
        Width of the sensitivity curve
    label : `string`
        Label for the sensitivity curve in legends
    **kwargs : `various`
        Keyword args are passed to :meth:`legwork.psd.power_spectral_density`, see those docs for details on
        possible arguments.
    Returns
    -------
    fig : `matplotlib Figure`
        The figure on which the distribution is plotted
    ax : `matplotlib Axis`
        The axis on which the distribution is plotted
    """
    if frequency_range is None:
        frequency_range = np.logspace(-5, 0, 1000) * u.Hz

    if fig is None or ax is None:
        fig, ax = plt.subplots()

    # work out what the noise amplitude should be
    PSD = psd.power_spectral_density(f=frequency_range, **kwargs)
    if y_quantity == "ASD":
        noise_amplitude = np.sqrt(PSD)
    elif y_quantity == "h_c":
        noise_amplitude = np.sqrt(frequency_range * PSD)
    else:
        raise ValueError("y_quantity must be one of 'ASD' or 'h_c'")

    # plot the curve and fill if needed
    with quantity_support():
        ax.loglog(frequency_range, noise_amplitude, color=color, label=label, linewidth=linewidth)
        if fill:
            ax.fill_between(frequency_range, np.zeros_like(noise_amplitude), noise_amplitude,
                            alpha=alpha, color=color)

    # adjust labels, sizes and frequency limits to plot is flush to the edges
    ax.set_xlabel(r'Frequency [$\rm Hz$]')
    if y_quantity == "ASD":
        ax.set_ylabel(r'ASD $[\rm Hz^{-1/2}]$')
    else:
        ax.set_ylabel(r'Characteristic Strain')

    ax.tick_params(axis='both', which='major')
    ax.set_xlim(np.min(frequency_range).value, np.max(frequency_range).value)

    if show:
        plt.show()

    return fig, ax