import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import astropy.units as u
import GravitationalWaves.psd as psd
from astropy.visualization import quantity_support

plt.rc('font', family='serif') # set the default font and fontsize
plt.rcParams['text.usetex'] = False
fs = 24
params = {'figure.figsize': (12, 8), 'legend.fontsize': fs, 'axes.labelsize': fs, 'xtick.labelsize': 0.7 * fs, 'ytick.labelsize': 0.7 * fs}
plt.rcParams.update(params) # update various fontsizes to match

__all__ = ['plot_1D_dist', 'plot_2D_dist', 'plot_sensitivity_curve']


def plot_1D_dist(x, weights=None, disttype="hist", fig=None, ax=None, xlabel=None, ylabel=None,
                 xlim=None, ylim=None, color=None, show=True, **kwargs):
    
    """
    plot a 1D distribution of ``x``.
    
    ### Parameters
        x : `float/int array` Variable to plot, should be a 1D array
        weights : `float/int array` Weights for each variable in ``x``, must have the same shape
        disttype : `{{ "hist", "kde", "ecdf" }}` Which type of distribution plot to use
        fig: `matplotlib Figure` A figure on which to plot the distribution. Both `ax` and `fig` must be supplied for either to be used
        ax: `matplotlib Axis` An axis on which to plot the distribution. Both `ax` and `fig` must be supplied for either to be used
        xlabel : `string` Label for the x axis, passed to Axes.set_xlabel()
        ylabel : `string` Label for the y axis, passed to Axes.set_ylabel()
        xlim : `tuple` Lower and upper limits for the x axis, passed to Axes.set_xlim()
        ylim : `tuple` Lower and upper limits for the y axis, passed to Axes.set_ylim()
        color : `string or tuple` Colour to use for the plot
        show : `boolean` Whether to immediately show the plot or only return the Figure and Axis
        **kwargs : `(if disttype=="hist")` Include values for any `bins, range, density, log, label` or more.
        **kwargs : `(if disttype=="kde")` Include values for any `gridsize, legend, fill, linewidth` or more.
        **kwargs : `(if disttype=="ecdf")` Include values for any `stat, linestyle, log_scale` or more.
    
    ### Returns:
        fig : `matplotlib Figure` The figure on which the distribution is plotted
        ax : `matplotlib Axis`The axis on which the distribution is plotted
    """
    
    if fig is None or ax is None: # create figure and axes
        fig, ax = plt.subplots()

    hist_args = {"bins": "auto", "density": True} # change default kwargs for matplotlib.hist
    kde_args = {"gridsize": 200, "legend": True} # change default kwargs for seaborn.kdeplot
    ecdf_args = {"stat": 'proportion', "legend": True} # change default kwargs for seaborn.ecdfplot

    plot_args = hist_args if disttype == "hist" else kde_args if disttype == "kde" else ecdf_args # set defaults
    
    for key, value in kwargs.items(): # update defaults
        plot_args[key] = value

    with quantity_support(): # create requested plot
        if disttype == "hist":
                ax.hist(x, weights=weights, color=color, **plot_args)
        elif disttype == "kde":
            sns.kdeplot(x=x, weights=weights, ax=ax, color=color, **plot_args)
        elif disttype == "ecdf":
            sns.ecdfplot(x=x, weights=weights, ax=ax, color=color, **plot_args)

    if xlabel is not None: # update axis labels
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
        
    if xlim is not None: # update axis limits
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if show: # show plot if requested
        plt.show()

    return fig, ax



def plot_2D_dist(x, y, weights=None, disttype="scatter", scatter_s=20, fig=None, ax=None,
                 xlabel=None, ylabel=None, xlim=None, ylim=None, color=None,
                 show=True, **kwargs):
    
    """
    Plot a 2D distribution of `x` and `y`
    
    ### Parameters
        x `float/int array`: Variable to plot on the x axis, should be a 1D array
        y : `float/int array` Variable to plot on the y axis, should be a 1D array
        weights : `float/int array` Weights for each variable pair (``x``, ``y``)
        disttype : `{{ "scatter", "kde" }}` Which type of distribution plot to use
        scatter_s : (`float`) Scatter point size, passed as ``s`` to a scatter plot and ignored for a KDE
        fig: `matplotlib Figure` A figure on which to plot the distribution. 
        ax: `matplotlib Axis` An axis on which to plot the distribution. 
        xlabel : `string` Label for the x axis
        ylabel : `string` Label for the y axis
        xlim : `tuple` Lower and upper limits for the x axis, passed to Axes.set_xlim()
        ylim : `tuple` Lower and upper limits for the y axis, passed to Axes.set_ylim()
        color : `string or tuple` Colour to use for the plot
        show : `boolean` Whether to immediately show the plot or only return the Figure and Axis
        **kwargs : `(if disttype=="hist")` Include values for any `bins, range, density, log, label` or more.
        **kwargs : `(if disttype=="kde")` Include values for any `gridsize, legend, fill, linewidth` or more.
        **kwargs : `(if disttype=="ecdf")` Include values for any `stat, linestyle, log_scale` or more.
    
    ### Returns:
        fig : `matplotlib Figure` The figure on which the distribution is plotted
        ax : `matplotlib Axis`The axis on which the distribution is plotted
    """
    
    if fig is None or ax is None: # create figure and axes
        fig, ax = plt.subplots()

    scatter_args = {} # change default kwargs for matplotlib.scatter
    kde_args = {"gridsize": 200, "legend": True, "levels": 10, "thresh": 0.02} # kwargs for seaborn.kdeplot
    plot_args = scatter_args if disttype == "scatter" else kde_args

    for key, value in kwargs.items(): # update values supplied
        plot_args[key] = value

    with quantity_support(): # create requested plot
        if disttype == "scatter":
            if weights is not None: # change size of points based on weights
                scatter_s = weights * scatter_s
            ax.scatter(x, y, s=scatter_s, color=color, **plot_args)
        elif disttype == "kde":
            sns.kdeplot(x=x, y=y, weights=weights, ax=ax, color=color, **plot_args)

    if xlabel is not None: # update axis labels
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
        
    if xlim is not None: # update axis limits
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if show: # show plot if requested
        plt.show()

    return fig, ax



def plot_sensitivity_curve(frequency_range=None, y_quantity="ASD", fig=None, ax=None, show=True,
                           color="#18068b", fill=True, alpha=0.2, linewidth=1, label=None, **kwargs):
    
    """
    Plot the LISA sensitivity curve
    
    ### Parameters
        frequency_range : `float array` Frequency values at which to plot the sensitivity curve
        y_quantity : `{{ "ASD", "h_c" }}` Which quantity to plot on y axis (amplitude spectral density or characteristic strain)
        fig: `matplotlib Figure` A figure on which to plot the distribution. 
        ax: `matplotlib Axis` An axis on which to plot the distribution. 
        show : `boolean` Whether to immediately show the plot or only return the Figure and Axis
        color : `string or tuple` Colour to use for the curve, 
        fill : `boolean` Whether to fill the area below the sensitivity curve
        alpha : `float` Opacity of the filled area below the sensitivity curve (ignored if fill is `False`)
        linewidth : `float` Width of the sensitivity curve
        label : `string` Label for the sensitivity curve in legends
        **kwargs : `various` Keyword args are passed
        
    ### Returns
        fig : `matplotlib Figure` The figure on which the distribution is plotted
        ax : `matplotlib Axis` The axis on which the distribution is plotted
    """
    
    if frequency_range is None:
        frequency_range = np.logspace(-5, 0, 1000) * u.Hz
    if fig is None or ax is None:
        fig, ax = plt.subplots()

    PSD = psd.power_spectral_density(f=frequency_range, **kwargs) # get noise amplitude
    if y_quantity == "ASD":
        noise_amplitude = np.sqrt(PSD)
    elif y_quantity == "h_c":
        noise_amplitude = np.sqrt(frequency_range * PSD)
    else:
        raise ValueError("y_quantity must be one of 'ASD' or 'h_c'")

    with quantity_support(): # plot curve and fill if needed
        ax.loglog(frequency_range, noise_amplitude, color=color, label=label, linewidth=linewidth)
        if fill:
            ax.fill_between(frequency_range, np.zeros_like(noise_amplitude), noise_amplitude, alpha=alpha, color=color)

    ax.set_xlabel(r'Frequency [$\rm Hz$]') # adjust labels, sizes and frequency limits
    if y_quantity == "ASD":
        ax.set_ylabel(r'ASD $[\rm Hz^{-1/2}]$')
    else:
        ax.set_ylabel(r'Characteristic Strain')

    ax.tick_params(axis='both', which='major')
    ax.set_xlim(np.min(frequency_range).value, np.max(frequency_range).value)

    if show:
        plt.show()

    return fig, ax