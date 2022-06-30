
import os
import numpy as np
import matplotlib.pyplot as plt
from pypeit.core import extract
from pypeit.core.flux_calib import Flam_to_Nlam
from pypeit.core.wavecal import wvutils
from astropy.io import ascii, fits
from scipy import interpolate
from astropy.stats import sigma_clipped_stats


from IPython import embed

from pypeit import slittrace
from pypeit import specobjs
from pypeit import io

from pypeit.display import display
from pypeit.images.imagebitmask import ImageBitMask
from pypeit.images.detector_container import DetectorContainer
from pypeit import masterframe
from pypeit import spec2dobj
from pypeit.scripts import scriptbase


from IPython import embed

def simulate_qso_Npix(redshift, exptime, debug=False):
    # TODO: use pypeit functions instead
    # extinction? telluric?

    # read in quasar spectrum from Selsing 2015, in erg/s/cm2/A
    # shift the spectrum to the observe frame and cut off at the Lyman-alpha
    # then convert to N_pixel
    
    dat = ascii.read("../resource/Selsing2015.dat")
    wl_rest = dat["col1"]
    wl_obs = wl_rest * (1 + redshift)
    flux = dat["col2"] # in 1e-17 erg/s/cm2/A
    flux_err = dat["col3"]

    # absorption of HI
    # TODO: simulate "real" absorption
    # wl_lya = 1215.67 * (1 + redshift)
    wl_lya = 0
    trough = wl_obs < wl_lya
    flux[trough] = 0

    # use one archived sensitivity function
    sens = fits.getdata("../resource/GD153_lris_sens.fits", ext=2)
    wl_sens = sens["SENS_WAVE"]
    wl_gpm = (wl_sens>1)
    zp = sens["SENS_ZEROPOINT"]
    zp_gpm = sens["SENS_ZEROPOINT_GPM"]
    zp_fit = sens["SENS_ZEROPOINT_FIT"]

    # calculate the 1/S_lambda, the factor that converts the flux to N_lambda
    factor = Flam_to_Nlam(wl_sens[wl_gpm], zp[wl_gpm]).flatten()
    func_factor = interpolate.interp1d(wl_sens[wl_gpm], factor)

    mask = (wl_obs > wl_sens[wl_gpm][0]) & (wl_obs < wl_sens[wl_gpm][-1])
    Nlam = flux[mask] * func_factor(wl_obs[mask])
    # Nlam_err = flux_err[mask] * func_factor(wl_obs[mask])
    delta_wave = wvutils.get_delta_wave(wl_obs[mask], (wl_obs[mask] > 1.0)).flatten()
    Npix = Nlam * exptime * delta_wave
    # Npix_err = Nlam_err * exptime * delta_wave

    if debug:
        fig, ax = plt.subplots(figsize=(12,6))
        ax.plot(wl_obs[mask], Npix)
        ax.set_xlim(wl_sens[wl_gpm][0], wl_sens[wl_gpm][-1])
        ax.set_xlabel(r"Wavelength [$\AA$]")
        ax.set_ylabel(r"$\rm N_{pix}$ [$\mathrm{photons/pixel}$]")
        plt.show()

    return wl_obs[mask], Npix

def plot_qso(redshift):
    dat = ascii.read("../resource/Selsing2015.dat")
    wl_rest = dat["col1"]
    wl_obs = wl_rest * (1 + redshift)
    flux = dat["col2"] # in 1e-17 erg/s/cm2/A
    flux_err = dat["col3"]
    wl_lya = 1215.67 * (1 + redshift)
    trough = wl_obs < wl_lya
    flux[trough] = 0

    fig, ax = plt.subplots(figsize=(12,6))
    ax.plot(wl_obs, flux, label="Selsing 2015", color="black")
    ax.fill_between(wl_obs, flux-flux_err, flux+flux_err, color="black", alpha=0.2)
    ax.set_xlim(7300, 10500)
    # ymin, ymax = ax.get_ylim()
    # ax.vlines(wl_lya, ymin, ymax)
    plt.show()
    return

#if __name__ == "__main__":
wave, Ncounts = simulate_qso_Npix(6.5, 20)

sci_path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/keck_lris_red_mark4/long_600_10000_d680/Science/'
spec2dfile = os.path.join(sci_path, 'spec2d_r220127_00123-J1209+0135_OFF_LRISr_20220127T145655.277.fits')
spec1dfile = spec2dfile.replace('spec2d', 'spec1d')
#spec1dfile = os.path.join(sci_path, 'spec1d_r220127_00123-J1209+0135_OFF_LRISr_20220127T145655.277.fits')

# Load it up -- NOTE WE ALLOW *OLD* VERSIONS TO GO FORTH
det = 1
detname = DetectorContainer.get_name(det)
spec2DObj = spec2dobj.Spec2DObj.from_file(spec2dfile, detname, chk_version=False)
sobjs = specobjs.SpecObjs.from_fitsfile(spec1dfile, chk_version=False)

offset = 100
sobjs_fake = sobjs[0].copy()
sobjs_fake.TRACE_SPAT = sobjs_fake.TRACE_SPAT - offset

gpm = spec2DObj.bpmmask == 0
sciimg = spec2DObj.sciimg

# Boxcar extract a new object at this location to get the boxcar wavelengths
# TODO think about base_var, count_scale and noise_floor
extract.extract_boxcar(spec2DObj.sciimg, spec2DObj.ivarmodel, gpm, spec2DObj.waveimg, spec2DObj.skymodel, sobjs_fake,
                       base_var=None, count_scale=None, noise_floor=None)
# Make a plot of the location where you extracted
image = (spec2DObj.sciimg - spec2DObj.skymodel) * gpm
mean, med, sigma = sigma_clipped_stats(image[spec2DObj.bpmmask == 0], sigma_lower=5.0,
                                       sigma_upper=5.0)
cut_min = mean - 1.0 * sigma
cut_max = mean + 4.0 * sigma
chname_skysub = 'skysub'
viewer, ch_skysub = display.show_image(image, chname=chname_skysub,
                                       waveimg=spec2DObj.waveimg,
                                       cuts=(cut_min, cut_max))
display.show_trace(viewer, ch_skysub, sobjs_fake.TRACE_SPAT, 'fake object', color='orange')

# 1) Interpolate the counts per second onto the sobjs_fake.BOX_WAVE

# 2) Create a 2d image with an object profile comprising a Gaussian centered on TRACE_SPAT

# 3) Force the extract_boxcar of the 2d image to equal to the total counts/per second/per real data pixel

# 4) As a simple test you can simply add the 2d image to the noisy 2d skysubtracted image (think we showed in ginga above)
# which will already give you an idea of what a high-z quasar looks like





# plot_qso(6)