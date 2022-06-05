import imp
from .plot_spec1d import plot_spec1d, plot_single, plot_series
from .plot_acq import plot_acq, plot_hist, plot_acq_and_hist
from .plot_extinction import plot_extinction

__all__ = ['plot_spec1d', 'plot_single', 'plot_series', 
           'plot_acq', 'plot_hist', 'plot_acq_and_hist',
           'plot_extinction']

import matplotlib.pyplot as plt
plt.style.use('science')