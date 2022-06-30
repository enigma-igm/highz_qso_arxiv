from gettext import find
import imp
import os
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np

from IPython import embed

__all__ = ["find_peak", "r_theta", "get_mosfire_acq", "get_mosfire_acq_proc", 
           "circle_mask",
           "naive_bkg_subtract", "sigma_clipping_bkg_subtract"]

def find_peak(image, sig_thresh=3.0, sig_clip=3.0, box_size=5, max_only=True):
    # TODO: spike bug
    # perhaps we can convolve first
    from photutils.detection import find_peaks
    mean, median, std = sigma_clipped_stats(image, sigma=sig_clip)
    threshold = median + (sig_thresh * std)
    tbl = find_peaks(image, threshold, box_size=box_size)
    if max_only:
        peak = tbl[tbl["peak_value"] == np.max(tbl["peak_value"])]
        peak = np.transpose((peak['x_peak'], peak['y_peak']))
        return peak[0]
    else:
        return np.transpose((tbl['x_peak'], tbl['y_peak']))

def r_theta(im, xc, yc):
    """r_theta - make a radius mask and return the radius rr and the angle phi for point (xc,yc)"""
    ny, nx = im.shape
    yp, xp = np.mgrid[0:ny,0:nx]
    yp = yp - yc
    xp = xp - xc
    rr = np.sqrt(np.power(yp,2.) + np.power(xp,2.))
    phi = np.arctan2(yp, xp)
    return(rr, phi)

def circle_mask(image, xc, yc, radius):
    rr, phi = r_theta(image, xc, yc)
    mask = (rr <= radius)
    return mask

def get_mosfire_acq(path, skyfile, objfile, cut_min=None, cut_max=None, sig_min=3.0, sig_max=3.0):
    sky_file = os.path.join(path, skyfile)
    obj_file = os.path.join(path, objfile)

    obj = fits.getdata(obj_file)
    sky = fits.getdata(sky_file)
    diff = obj - sky

    # Get the sky level in the window
    ny, nx = diff.shape
    # define crude y, x boundaries
    x_coord, y_coord = np.meshgrid(np.arange(nx), np.arange(ny))

    # hardcode for mosfire acquisition box
    upper_left_y, upper_left_x = (1070, 1034)
    upper_right_y, upper_right_x = (1070, 1051)
    lower_left_y, lower_left_x = (1035, 1037)
    lower_right_y, lower_right_x = (1032, 1055)

    median_box = (x_coord > lower_left_x) & (x_coord < upper_right_x) & (y_coord > lower_left_y) & (y_coord < upper_right_y)
    mean_sky, med_sky, sigma_sky = sigma_clipped_stats(diff[median_box], sigma_lower=3.0, sigma_upper=3.0)
    cut_min = med_sky - sig_min*sigma_sky if cut_min is None else cut_min
    cut_max = med_sky + sig_max*sigma_sky if cut_max is None else cut_max

    image = diff * median_box

    image_cutout = image[1035+1:1070,1037+1:1051]
    return image_cutout

def get_mosfire_acq_proc(path, skyfile, objfile):
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit.images import buildimage
    sky_file = os.path.join(path, skyfile)
    obj_file = os.path.join(path, objfile)
    spectrograph = load_spectrograph("keck_mosfire")
    par = spectrograph.default_pypeit_par()['calibrations']['biasframe']
    det = 1
    img = buildimage.buildimage_fromlist(spectrograph, det, par,
                                         [obj_file], mosaic=False)
    bkg_img = buildimage.buildimage_fromlist(spectrograph, det, par, 
                                            [sky_file], mosaic=False)
    img = img.sub(bkg_img, par['process'])
    import matplotlib.pyplot as plt
    embed()
    # hardcode acq box here
    return img.image[1031:1050,1029:1063], img.ivar[1031:1050,1029:1063]

def naive_bkg_subtract(image, method="median", mask_source=False, mask_radius=6):
    """Naive background subtraction
       Subtracting a scalar, either median or biweight of the image

    Args:
        image (_type_): _description_
        method (str, optional): _description_. Defaults to "median".

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    # TODO: uncertainty
    if mask_source:
        peak = find_peak(image)
        mask = circle_mask(image, peak[0], peak[1], mask_radius)
    else:
        # array with all elements == False
        mask = (image == True)
    if method == "median":
        # return np.median(image)
        return image - np.median(image*~mask)
    elif method == "biweight":
        from astropy.stats import biweight_location
        # return biweight_location(image)
        return image - biweight_location(image*~mask)
    else:
        raise ValueError("No such method!")

def sigma_clipping_bkg_subtract(image, sigma=3.0, mask_source=True, mask_radius=5, show_mask=False):
    if mask_source:
        # from photutils.segmentation import make_source_mask
        # mask = make_source_mask(image, nsigma=2, npixels=5, dilate_size=11)
        peak = find_peak(image)
        mask = ~circle_mask(image, peak[0], peak[1], mask_radius)
        if show_mask:
            from ..plot import plot_acq
            plot_acq(image*mask, display=True)
    else:
        mask = None
    mean, median, std = sigma_clipped_stats(image, sigma=sigma, mask=mask)
    # return mean, median
    return image - median

