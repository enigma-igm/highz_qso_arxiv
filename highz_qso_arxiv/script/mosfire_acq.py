"""
Show difference MOSFIRE acquisition frames.
"""

# modified from: https://github.com/pypeit/PypeIt-development-suite/blob/master/dev_algorithms/mosfire/mosfire_acq.py

import sys
import os
import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import argparse

from pyds9 import DS9
from IPython import embed

def get_acq_image(path, prefix, obj_frame, sky_frame, cut_min=None, cut_max=None, sig_min=3.0, sig_max=3.0, plot_ds9=True):
    sky_file = os.path.join(path, '{:s}{:04d}.fits'.format(prefix, sky_frame))
    obj_file = os.path.join(path, '{:s}{:04d}.fits'.format(prefix, obj_frame))

    obj = fits.getdata(obj_file)
    sky = fits.getdata(sky_file)
    diff = obj - sky
    
    # Get the sky level in the window
    ny, nx = diff.shape
    # define crude y, x boundaries
    x_coord, y_coord = np.meshgrid(np.arange(nx), np.arange(ny))
    upper_left_y, upper_left_x = (1070, 1034)
    upper_right_y, upper_right_x = (1070, 1051)
    lower_left_y, lower_left_x = (1035, 1037)
    lower_right_y, lower_right_x = (1032, 1055)

    median_box = (x_coord > lower_left_x) & (x_coord < upper_right_x) & (y_coord > lower_left_y) & (y_coord < upper_right_y)
    mean_sky, med_sky, sigma_sky = sigma_clipped_stats(diff[median_box], sigma_lower=3.0, sigma_upper=3.0)
    cut_min = med_sky - sig_min*sigma_sky if cut_min is None else cut_min
    cut_max = med_sky + sig_max*sigma_sky if cut_max is None else cut_max

    image = diff * median_box
    if plot_ds9:
        try:
            d = DS9()
        except TypeError:
            print("Perhaps you forget to turn DS9 on!")
            exit()
        d.set_np2arr(image)
        d.set('scale limits {:5.3} {:5.3}'.format(cut_min, cut_max))
        d.set('scale linear')
        d.set('zoom to fit')
        d.set('zoom 5')
    return image

def main():
    path = "../resource"
    get_acq_image(path, "qso_", 1, 0)

if __name__ == '__main__':
    main()