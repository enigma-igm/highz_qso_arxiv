
import os
import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from scipy import interpolate
from astropy.stats import sigma_clipped_stats
from astropy.cosmology import FlatLambdaCDM, z_at_value

from highz_qso_arxiv.resource.filters import ukirt
from highz_qso_arxiv.util.spec2dutil import gauss_comb
from highz_qso_arxiv.util import luminosity_to_flux, redshift_to_distance
from highz_qso_arxiv.util.spec1dutil import gp_trough, damping_wing, extend_to_lower

from pypeit import io
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit.core import extract
from pypeit.display import display
from pypeit.sensfunc import IRSensFunc
from pypeit.core.wavecal import wvutils
from pypeit.core.moment import moment1d
from pypeit.core.flux_calib import Flam_to_Nlam
from pypeit.images.detector_container import DetectorContainer

from IPython import embed

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def load_quasar():
    # TODO: unit of wavelength
    dat = ascii.read("../resource/Selsing2015.dat")
    wl_rest, flux = dat["col1"], dat["col2"] * 1e-17 * u.erg / u.s / u.cm**2 / u.AA
    return wl_rest, flux

def load_star():
    dat = ascii.read(f"../resource/dwarf/L1_2MASS0746+20.dat")
    wl_obs, flux = dat["col1"], dat["col2"] * u.erg / u.s / u.cm**2 / u.AA
    return wl_obs, flux

def parse_quasar(wave_rest, flux, redshift, m_1450=None, M_1450=None, m_J=None):
    """
    Parse quasar spectrum
    """
    dL = redshift_to_distance(redshift, cosmo)
    wave_obs = wave_rest * (1 + redshift)

    # for scaling the qso template
    mask_1450 = (wave_rest > 1445) & (wave_rest < 1455)
    flam_1450_temp = np.median(flux[mask_1450])

    if m_1450 != None:
        mAB = m_1450 * u.ABmag
        flam_1450 = mAB.to(1e-17 * u.erg/u.s/u.cm**2/u.AA, u.spectral_density(1450*u.AA*(1+redshift)))
        scale = flam_1450 / flam_1450_temp
    elif M_1450 != None:
        mAB = (M_1450 + 5 * np.log10(dL*1e6/10)) * u.ABmag
        flam_1450 = mAB.to(1e-17 * u.erg/u.s/u.cm**2/u.AA, u.spectral_density(1450*u.AA*(1+redshift)))
        scale = flam_1450 / flam_1450_temp
    elif m_J != None:
        m_J_temp = ukirt.get_ab_magnitudes(flux, wave_obs*u.AA)[0][0]
        scale = 10**(-(m_J-m_J_temp)/2.5)
    else:
        raise ValueError('m_1450 or M_1450 or m_J must be specified')
    flux = flux * scale
    wave_obs, flux = extend_to_lower(wave_obs, flux, 5000.)
    flux = gp_trough(wave_obs, flux, redshift)
    return wave_obs, flux

def parse_star(wave_rest, flux, m_J):
    """
    Parse star spectrum
    """
    m_J_temp = ukirt.get_ab_magnitudes(flux, wave_rest*u.AA)[0][0]
    scale = 10**(-(m_J-m_J_temp)/2.5)
    flux = flux * scale
    return wave_rest, flux

def flux_to_Nlam(flux, wave_obs, sens):
    """
    Convert flux to Nlam
    """
    # calculate the 1/S_lambda, the factor that converts the flux to N_lambda
    S_lam_units = 1e-17 * u.erg / u.cm**2
    factor = Flam_to_Nlam(sens.wave, sens.zeropoint).flatten() 
    facotr_interp = interpolate.interp1d(sens.wave.flatten(), factor)

    mask = (wave_obs > sens.wave[0]) & (wave_obs < sens.wave[-1])
    Nlam = flux[mask] * facotr_interp(wave_obs[mask]) / S_lam_units
    Nlam = Nlam.to(1/u.AA/u.s).value # 1 / (Angstrom s)
    return wave_obs[mask], Nlam

def simulate(sens, spec2DObj, sobjs, trace_id, offset, exptime, load_func, parse_func, show_trace=False, **kwargs):
    wl, flux = load_func()
    wl, flux = parse_func(wl, flux, **kwargs)
    wl, Nlam = flux_to_Nlam(flux, wl, sens)

    sobjs_fake = sobjs[trace_id].copy()
    sobjs_fake.TRACE_SPAT = sobjs_fake.TRACE_SPAT + offset

    gpm = spec2DObj.bpmmask == 0

    # Boxcar extract a new object at this location to get the boxcar wavelengths
    # TODO think about base_var, count_scale and noise_floor
    extract.extract_boxcar(spec2DObj.sciimg, spec2DObj.ivarmodel, gpm, spec2DObj.waveimg, spec2DObj.skymodel, sobjs_fake,
                        base_var=None, count_scale=None, noise_floor=None)

    # Interpolate Nlam onto the sobjs_fake.BOX_WAVE
    Nlam_interp = interpolate.interp1d(wl, Nlam)
    mask = (sobjs_fake.BOX_WAVE > wl[0]) & (sobjs_fake.BOX_WAVE < wl[-1])
    try:
        Nlam_boxwave = Nlam_interp(sobjs_fake.BOX_WAVE)
    except ValueError:
        raise ValueError("Nlam_interp failed")
    delta_wave = wvutils.get_delta_wave(sobjs_fake.BOX_WAVE, (sobjs_fake.BOX_WAVE > 1.0))
    Npix_boxwave = Nlam_boxwave * exptime * delta_wave

    # Create a 2d image with an object profile comprising a Gaussian centered on TRACE_SPAT
    img_gaussians = gauss_comb(spec2DObj.sciimg.shape, sobjs_fake.TRACE_SPAT, sig=sobjs_fake.FWHM/2.355)

    # Force the extract_boxcar of the 2d image to equal to the total counts/per second/per real data pixel
    flux_gaussians = moment1d(img_gaussians, sobjs_fake.TRACE_SPAT, sobjs_fake.BOX_RADIUS * 2, order=[0])[0]
    factor_fake = Npix_boxwave / flux_gaussians
    img_fake = (img_gaussians.T * factor_fake).T
    print('FWHM:', sobjs_fake.FWHM)

    # Make a plot of fake quasar
    image = (spec2DObj.sciimg - spec2DObj.skymodel) * gpm
    mean, med, sigma = sigma_clipped_stats(image[spec2DObj.bpmmask == 0], 
                                           sigma_lower=5.0,
                                           sigma_upper=5.0)
    cut_min, cut_max = mean - 1.0 * sigma, mean + 4.0 * sigma
    viewer, ch_skysub = display.show_image(image+img_fake, chname='skysub',
                                           waveimg=spec2DObj.waveimg,
                                           cuts=(cut_min, cut_max))
    if show_trace:
        display.show_trace(viewer, ch_skysub, sobjs_fake.TRACE_SPAT, 'fake object', color='orange')
    return img_fake

# hdul = io.fits_open("../resource/sensfunc/GD153_lris_long_8.7_sens.fits")
# hdul = io.fits_open("../resource/sensfunc/GD153_mosfire_sens.fits")
hdul = io.fits_open("../resource/sensfunc/GD153_lris_sens.fits")
sens = IRSensFunc.from_hdu(hdul)

sci_path = '../arxiv/LRIS_2203/LRIS_220306/reduced/sim/Science/'
# sci_path = '../arxiv/MOSFIRE_2201/reduced/all/Science'

spec2dfile = os.path.join(sci_path, 'spec2d_r220306_00045-J0901+2906_OFF_LRISr_20220306T055449.306.fits')
# spec2dfile = os.path.join(sci_path, 'spec2d_MF.20220111.33636-J0637+3812_MOSFIRE_20220111T092036.865.fits')
spec1dfile = spec2dfile.replace('spec2d', 'spec1d')

detname = DetectorContainer.get_name(det=1)
spec2DObj = spec2dobj.Spec2DObj.from_file(spec2dfile, detname, chk_version=False)
sobjs = specobjs.SpecObjs.from_fitsfile(spec1dfile, chk_version=False)

simulate(sens=sens, spec2DObj=spec2DObj, sobjs=sobjs, trace_id=2, offset=-100, exptime=300., 
         load_func=load_quasar, parse_func=parse_quasar, show_trace=False, redshift=7.5, m_1450=21.5)
# simulate(sens=sens, spec2DObj=spec2DObj, sobjs=sobjs, trace_id=2, offset=-100, exptime=300., 
#          load_func=load_star, parse_func=parse_star, show_trace=False, m_J=21.3)