import imp
from .plot_spec1d import plot_spec1d, plot_single, plot_series, plot_template_dat
from .plot_acq import plot_acq, plot_hist, plot_acq_and_hist
from .plot_extinction import plot_extinction
from .plot_contour import plot_contour
from .plot_covariance import plot_cov, error_cov
from .plot_colorline import plot_cline

__all__ = ['plot_spec1d', 'plot_single', 'plot_series', 
           'plot_acq', 'plot_hist', 'plot_acq_and_hist',
           'plot_extinction', 'plot_template_dat', 'plot_contour', 'plot_cov', 'error_cov', 'plot_cline']

import matplotlib.pyplot as plt
plt.style.use('science')

import matplotlib as mpl
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 