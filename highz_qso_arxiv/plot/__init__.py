import imp
from .plot_spec1d import plot_spec1d, plot_single, plot_series, plot_template_dat
from .plot_acq import plot_acq, plot_hist, plot_acq_and_hist
from .plot_extinction import plot_extinction
from .plot_contour import plot_contour

__all__ = ['plot_spec1d', 'plot_single', 'plot_series', 
           'plot_acq', 'plot_hist', 'plot_acq_and_hist',
           'plot_extinction', 'plot_template_dat', 'plot_contour']

import matplotlib.pyplot as plt
plt.style.use('science')