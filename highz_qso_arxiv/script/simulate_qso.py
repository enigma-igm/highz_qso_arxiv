
import os
from re import U
import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from pypeit.core import extract
from pypeit.core.wavecal import wvutils
from pypeit.core.moment import moment1d
from pypeit.core.flux_calib import Flam_to_Nlam
from pypeit.sensfunc import IRSensFunc
from astropy.io import ascii, fits
from scipy import interpolate
from astropy.stats import sigma_clipped_stats
from astropy.cosmology import FlatLambdaCDM, z_at_value

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
from scipy.special import erf

from IPython import embed

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def gauss_comb(shape, center, sig=5.):
    img = np.zeros(shape, dtype=float)
    x = np.arange(shape[1])
    for i,_c in enumerate(center):
        img[i,:] = (erf((x-center[i]+0.5)/np.sqrt(2)/sig) 
                   - erf((x-center[i]-0.5)/np.sqrt(2)/sig))/2.
    return img

def zfill(wave, counts, wave_min):
    """
    expand wave to wave_min, fill in zero
    """
    if wave[0] < wave_min:
        return wave, counts
    wave_min_old = np.min(wave)
    wave_new = np.arange(wave_min, wave_min_old, 1)
    counts_new = np.zeros_like(wave_new)
    wave_new = np.append(wave_new, wave)
    counts_new = np.append(counts_new, counts)
    return wave_new, counts_new

def luminosity_to_flux(luminosity, luminosity_distance):
    """
    convert luminosity to flux
    """
    return luminosity / (4 * np.pi * luminosity_distance**2)

def redshift_to_distance(redshift) -> float:
    """
    Returns the distance in Mpc for a given redshift.
    """
    return cosmo.luminosity_distance(redshift).value
    
def simulate_Nlam(model='qso', redshift=None, M_1450=None, 
                  star_scale=None, star_type=None, debug=False):
    # read in quasar spectrum from Selsing 2015, in erg/s/cm2/A
    # shift the spectrum to the observe frame and cut off at the Lyman-alpha
    # then convert to Nlam
    
    # TODO 
    # 1. use speclite to resample
    # 2. scale for star template (use magnitude)
    # 3. quasar absorption

    if model == 'qso':        
        dL = redshift_to_distance(redshift)
        mAB = M_1450 + 5 * np.log10(dL*1e6/10)
        mAB = mAB * u.ABmag
        flam = mAB.to(u.erg/u.s/u.cm**2/u.AA, u.spectral_density(1450*u.AA))

        dat = ascii.read("../resource/Selsing2015.dat")
        wl_rest = dat["col1"]
        wl_obs = wl_rest * (1 + redshift)
        flux = dat["col2"] # in 1e-17 erg/s/cm2/A
        flux_err = dat["col3"]

        mask_1450 = (wl_rest > 1445) & (wl_rest < 1455)
        flux_mean = np.mean(flux[mask_1450])
        scale = flam / flux_mean
        flux = (flux * scale * 1e17).value
        # absorption of HI
        # TODO: simulate "real" absorption
        wl_lya = 1215.67 * (1 + redshift)
        # wl_lya = 0
        trough = wl_obs < wl_lya
        flux[trough] = 0
        
    elif model == 'star':
        dat = ascii.read(f"../resource/dwarf/keck_lris_{star_type}_1.dat")
        wl_obs = dat["col1"]
        flux = dat["col2"] * 1e17 * star_scale # in 1e-17 erg/s/cm2/A

    # use one archived sensitivity function
    hdul = io.fits_open("../resource/GD153_lris_sens.fits")
    sens = IRSensFunc.from_hdu(hdul)
    # wl_gpm = (sens.wave>1)

    # calculate the 1/S_lambda, the factor that converts the flux to N_lambda
    factor = Flam_to_Nlam(sens.wave, sens.zeropoint).flatten()
    func_factor = interpolate.interp1d(sens.wave.flatten(), factor)

    mask = (wl_obs > sens.wave[0]) & (wl_obs < sens.wave[-1])
    wl_obs_lris = wl_obs[mask]
    Nlam = flux[mask] * func_factor(wl_obs_lris)
    
    # TODO this is too twisted, should generalize this
    wl_lris, Nlam_lris = zfill(wl_obs_lris, Nlam, sens.wave[0])

    return wl_lris, Nlam_lris

exptime = 20 # seconds
# wave, Nlam = simulate_Nlam(model='qso', redshift=6.1, M_1450=-26)
wave, Nlam = simulate_Nlam(model='star', star_scale=0.04, star_type='L0')

# sci_path = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/keck_lris_red_mark4/long_600_10000_d680/Science/'
sci_path = '../arxiv/LRIS_2203/LRIS_220306/reduced/all/Science/'
# spec2dfile = os.path.join(sci_path, 'spec2d_r220127_00123-J1209+0135_OFF_LRISr_20220127T145655.277.fits')
spec2dfile = os.path.join(sci_path, 'spec2d_r220306_00045-J0901+2906_OFF_LRISr_20220306T055449.306.fits')
spec1dfile = spec2dfile.replace('spec2d', 'spec1d')
#spec1dfile = os.path.join(sci_path, 'spec1d_r220127_00123-J1209+0135_OFF_LRISr_20220127T145655.277.fits')

# Load it up -- NOTE WE ALLOW *OLD* VERSIONS TO GO FORTH
det = 1
detname = DetectorContainer.get_name(det)
spec2DObj = spec2dobj.Spec2DObj.from_file(spec2dfile, detname, chk_version=False)
sobjs = specobjs.SpecObjs.from_fitsfile(spec1dfile, chk_version=False)

offset = 100
sobjs_fake = sobjs[1].copy()
sobjs_fake.TRACE_SPAT = sobjs_fake.TRACE_SPAT - offset
# sobjs_fake.BOX_WAVE vs. sobjs_fake.TRACE_SPAT

gpm = spec2DObj.bpmmask == 0
sciimg = spec2DObj.sciimg
# Boxcar extract a new object at this location to get the boxcar wavelengths
# TODO think about base_var, count_scale and noise_floor
extract.extract_boxcar(spec2DObj.sciimg, spec2DObj.ivarmodel, gpm, spec2DObj.waveimg, spec2DObj.skymodel, sobjs_fake,
                       base_var=None, count_scale=None, noise_floor=None)

# Make a plot of the location where you extracted
# image = (spec2DObj.sciimg - spec2DObj.skymodel) * gpm
# mean, med, sigma = sigma_clipped_stats(image[spec2DObj.bpmmask == 0], sigma_lower=5.0,
#                                        sigma_upper=5.0)
# cut_min = mean - 1.0 * sigma
# cut_max = mean + 4.0 * sigma
# chname_skysub = 'skysub'
# viewer, ch_skysub = display.show_image(image, chname=chname_skysub,
#                                        waveimg=spec2DObj.waveimg,
#                                        cuts=(cut_min, cut_max))
# display.show_trace(viewer, ch_skysub, sobjs_fake.TRACE_SPAT, 'fake object', color='orange')

# 1) Interpolate the counts per second onto the sobjs_fake.BOX_WAVE

# 2) Create a 2d image with an object profile comprising a Gaussian centered on TRACE_SPAT

# 3) Force the extract_boxcar of the 2d image to equal to the total counts/per second/per real data pixel

# 4) As a simple test you can simply add the 2d image to the noisy 2d skysubtracted image (think we showed in ginga above)
# which will already give you an idea of what a high-z quasar looks like

func_Nlamb = interpolate.interp1d(wave, Nlam)
mask = (sobjs_fake.BOX_WAVE > wave[0]) & (sobjs_fake.BOX_WAVE < wave[-1])
Nlam_boxwave = func_Nlamb(sobjs_fake.BOX_WAVE[mask])
delta_wave = wvutils.get_delta_wave(sobjs_fake.BOX_WAVE[mask], (sobjs_fake.BOX_WAVE[mask] > 1.0))
Npix_boxwave = Nlam_boxwave * exptime * delta_wave

img_gaussians = gauss_comb(spec2DObj.sciimg[mask].shape, sobjs_fake.TRACE_SPAT[mask], sig=1)
flux_gaussians = moment1d(img_gaussians, sobjs_fake.TRACE_SPAT[mask], sobjs_fake.BOX_RADIUS*2, order=[0])[0]
factor_fake = Npix_boxwave / flux_gaussians
img_fake_tmp = (img_gaussians.T*factor_fake).T

# Make a plot of fake quasar
image = (spec2DObj.sciimg - spec2DObj.skymodel) * gpm
img_fake = np.zeros_like(image)
img_fake[mask] = img_fake_tmp

mean, med, sigma = sigma_clipped_stats(image[spec2DObj.bpmmask == 0], sigma_lower=5.0,
                                       sigma_upper=5.0)
cut_min = mean - 1.0 * sigma
cut_max = mean + 4.0 * sigma
chname_skysub = 'skysub'
viewer, ch_skysub = display.show_image(image+img_fake, chname=chname_skysub,
                                       waveimg=spec2DObj.waveimg,
                                       cuts=(cut_min, cut_max))
# display.show_trace(viewer, ch_skysub, sobjs_fake.TRACE_SPAT, 'fake object', color='orange')