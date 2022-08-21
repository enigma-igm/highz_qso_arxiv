
from cProfile import label
import os
import speclite
import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from scipy import interpolate
from astropy.stats import sigma_clipped_stats
from astropy.cosmology import FlatLambdaCDM, z_at_value

from highz_qso_arxiv.resource.filters import ukirt_J, hsc_z, decam_z, sdss_i
from highz_qso_arxiv.util.spec2dutil import gauss_comb
from highz_qso_arxiv.util import luminosity_to_flux, redshift_to_distance, get_project_root, ivarsmooth
from highz_qso_arxiv.util.spec1dutil import add_gp_trough, add_damping_wing, add_telluric, extend_to_lower

from pypeit import io
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit.utils import inverse
from pypeit.core import extract
from pypeit.display import display
from pypeit.sensfunc import IRSensFunc
from pypeit.core.wavecal import wvutils
from pypeit.core.moment import moment1d
from pypeit.core.flux_calib import Flam_to_Nlam
from pypeit.spectrographs.util import load_spectrograph
from pypeit.images.detector_container import DetectorContainer
from pypeit.core import procimg

from IPython import embed

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

root = get_project_root()

def load_quasar():
    # TODO: unit of wavelength
    dat = ascii.read(os.path.join(root, 'resource/Selsing2015.dat'))
    wl_rest, flux = dat["col1"], dat["col2"] * 1e-17 * u.erg / u.s / u.cm**2 / u.AA
    return wl_rest, flux

def load_star():
    dat = ascii.read(os.path.join(root, f'resource/dwarf/combine_lris_M9.dat'))
    wl_obs, flux = dat["wave"], dat["flux"] * u.erg / u.s / u.cm**2 / u.AA
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
        m_J_temp = ukirt_J.get_ab_magnitudes(flux, wave_obs*u.AA)[0][0]
        scale = 10**(-(m_J-m_J_temp)/2.5)
    else:
        raise ValueError('m_1450 or M_1450 or m_J must be specified')
    flux = flux * scale
    wave_obs, flux = extend_to_lower(wave_obs, flux, 5000.)
    flux = add_gp_trough(wave_obs, flux, redshift)
    return wave_obs, flux

def parse_star(wave_rest, flux, m_J):
    """
    Parse star spectrum
    """
    # flux, wl = sdss_i.pad_spectrum(flux.value, wave_rest.value, method='zero')
    flux, wl = ukirt_J.pad_spectrum(flux.value, wave_rest.value, method='zero')
    # m_i = m_J + 1.99 + 0.82 + 0.96 + 0.37 - 0.938 # testing
    # m_i_temp = sdss_i.get_ab_magnitudes(flux*u.erg/u.s/u.cm**2/u.AA, wl*u.AA)[0][0]
    m_J_temp = ukirt_J.get_ab_magnitudes(flux, wl*u.AA)[0][0]
    scale = 10**(-(m_J-m_J_temp)/2.5)
    # scale = 10**(-(m_i-m_i_temp)/2.5)
    flux = flux * scale * u.erg/u.s/u.cm**2/u.AA
    return wl, flux

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

def simple_slit_loss(slitwidth, fwhm, platescale, binning):
    from scipy.special import erf
    return erf((slitwidth/2)/(np.sqrt(2)*fwhm*platescale*binning/2.355))

def simulate(sens, spec2DObj, sobjs, header, 
             trace_id, offset, exptime, load_func, parse_func, 
             telluric=None, damping_wing=None, slitloss=True,
             show_trace=False, debug=False, verbose=True, **kwargs):
    info = {}
    sobjs_fake = sobjs[trace_id].copy()
    sobjs_fake.TRACE_SPAT = sobjs_fake.TRACE_SPAT + offset
    sobjs_fake_noiseless = sobjs_fake.copy()


    wave, flux = load_func()
    wave, flux = parse_func(wave, flux, **kwargs)
    wl, Nlam = flux_to_Nlam(flux, wave, sens)
    if slitloss:
        slit_trans = simple_slit_loss(1, sobjs_fake.FWHMFIT, 0.123, 2)
        print('Median Slit transmission:', np.median(slit_trans))
        func_slit_trans = interpolate.interp1d(sobjs_fake.BOX_WAVE, slit_trans, fill_value="extrapolate")
        Nlam = Nlam * func_slit_trans(wl)
    if telluric is not None:
        wl, Nlam = add_telluric(wl, Nlam, telluric)

    gpm = spec2DObj.bpmmask == 0

    # Boxcar extract a new object at this location to get the boxcar wavelengths
    # TODO think about base_var, count_scale and noise_floor
    extract.extract_boxcar(spec2DObj.sciimg, spec2DObj.ivarmodel, gpm, spec2DObj.waveimg, spec2DObj.skymodel, sobjs_fake,
                           base_var=None, count_scale=None, noise_floor=None)

    # Interpolate Nlam onto the sobjs_fake.BOX_WAVE
    Nlam_interp = interpolate.interp1d(wl, Nlam, fill_value=0, bounds_error=False)
    mask = (sobjs_fake.BOX_WAVE > wl[0]) & (sobjs_fake.BOX_WAVE < wl[-1])
    try:
        Nlam_boxwave = Nlam_interp(sobjs_fake.BOX_WAVE)
    except ValueError:
        embed()
        raise ValueError("Nlam_interp failed")
    delta_wave = wvutils.get_delta_wave(sobjs_fake.BOX_WAVE, (sobjs_fake.BOX_WAVE > 1.0))
    Npix_boxwave = Nlam_boxwave * exptime * delta_wave

    # Create a 2d image with an object profile comprising a Gaussian centered on TRACE_SPAT
    img_gaussians = gauss_comb(spec2DObj.sciimg.shape, sobjs_fake.TRACE_SPAT, sig=sobjs_fake.FWHMFIT/2.355)

    # Force the extract_boxcar of the 2d image to equal to the total counts/per second/per real data pixel
    flux_gaussians = moment1d(img_gaussians, sobjs_fake.TRACE_SPAT, sobjs_fake.BOX_RADIUS * 2, order=[0])[0]
    factor_fake = Npix_boxwave / flux_gaussians
    img_fake = (img_gaussians.T * factor_fake).T
    if verbose:
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

    # before calculating signal-to-noise, we first define the region used in the calculation
    try:
        wave_min = (1 + kwargs['redshift']) * 1220
        wave_max = (1 + kwargs['redshift']) * 1230
    except KeyError:
        wave_min = 9500
        wave_max = 9600
    if verbose:
        print("wave_min:", wave_min)
        print("wave_max:", wave_max)
    snr_mask_fake = (sobjs_fake.BOX_WAVE < wave_max) & (sobjs_fake.BOX_WAVE > wave_min)
    snr_N_fake = len(sobjs_fake.BOX_WAVE[snr_mask_fake])
    snr_mask_target = (sobjs[trace_id].BOX_WAVE < wave_max) & (sobjs[trace_id].BOX_WAVE > wave_min)
    snr_N_target = len(sobjs[trace_id].BOX_WAVE[snr_mask_target])
    if verbose:
        print('SNR N:', snr_N_fake)

    # extract again to update sobjs_fake
    _base_var = inverse(spec2DObj.ivarmodel) - (spec2DObj.skymodel) - (spec2DObj.objmodel)
    _count_scale = None
    #adderr = 0.01
    #embed()
    #var = procimg.variance_model(_base_var, counts=spec2DObj.skymodel, count_scale=_count_scale, noise_floor=adderr)
    #ivarmodel = inverse(var)
    #rel_diff = np.mean((ivarmodel - spec2DObj.ivarmodel)/spec2DObj.ivarmodel)
    # JFH new line. Compare this to your noise realization
    extract.extract_boxcar(spec2DObj.skymodel+img_fake, spec2DObj.ivarmodel, gpm, spec2DObj.waveimg, spec2DObj.skymodel,
                           sobjs_fake_noiseless, base_var=_base_var, count_scale=None, noise_floor=None)
    extract.extract_boxcar(spec2DObj.sciimg+img_fake, spec2DObj.ivarmodel, gpm, spec2DObj.waveimg, spec2DObj.skymodel,
                           sobjs_fake, base_var=_base_var, count_scale=None, noise_floor=None)
    snr_signal_fake = np.sum(sobjs_fake_noiseless.BOX_COUNTS[snr_mask_fake]) / snr_N_fake
    snr_variance_fake = np.sum(inverse(sobjs_fake_noiseless.BOX_COUNTS_IVAR[snr_mask_fake])) / snr_N_fake**2
    snr_fake = snr_signal_fake / np.sqrt(snr_variance_fake)
    snr_signal_target = np.sum(sobjs[trace_id].BOX_COUNTS[snr_mask_target]) / snr_N_target
    snr_variance_target = np.sum(inverse(sobjs[trace_id].BOX_COUNTS_IVAR[snr_mask_target])) / snr_N_target**2
    snr_target = snr_signal_target / np.sqrt(snr_variance_target)
    if verbose:
        print('SNR_fake:', snr_fake)
        print('SNR_target:', snr_target)
    info['SNR'] = snr_fake

    if show_trace:
        display.show_trace(viewer, ch_skysub, sobjs_fake.TRACE_SPAT, 'fake object', color='orange')
    if debug:
        # sanity check - 1: noiseless spectrum vs. noisy spectrum
        plt.figure(figsize=(20,6))
        plt.plot(sobjs_fake.BOX_WAVE, sobjs_fake.BOX_COUNTS, drawstyle='steps-mid', color='black', label='noisy')
        plt.plot(sobjs_fake_noiseless.BOX_WAVE, sobjs_fake_noiseless.BOX_COUNTS, drawstyle='steps-mid', color='blue', label='noiseless')
        plt.plot(sobjs_fake.BOX_WAVE, np.sqrt(inverse(sobjs_fake.BOX_COUNTS_IVAR)), drawstyle='steps-mid', color='red')
        plt.plot(sobjs_fake.BOX_WAVE, sobjs_fake.BOX_COUNTS_SIG_DET, drawstyle='steps-mid', color='orange')
        plt.ylim(np.median(sobjs_fake.BOX_COUNTS)-3*np.median(np.sqrt(inverse(sobjs_fake.BOX_COUNTS_IVAR))), 
                 np.median(sobjs_fake.BOX_COUNTS)+8*np.median(np.sqrt(inverse(sobjs_fake.BOX_COUNTS_IVAR))))
        plt.title('noiseless vs. noisy spectrum', fontsize=40)
        plt.tick_params(labelsize=20)
        plt.xlabel('Wavelength [Angstrom]', fontsize=40)
        plt.ylabel('Counts', fontsize=40)
        plt.legend(fontsize=20)
        plt.show()

        # sanity check - 2: noise check
        from scipy.stats import norm

        plt.figure(figsize=(6,6))
        chi = (sobjs_fake.BOX_COUNTS-sobjs_fake_noiseless.BOX_COUNTS)/np.sqrt(inverse(sobjs_fake.BOX_COUNTS_IVAR))
        chi_mask = (chi < 10) & (chi > -10) 
        mean, std = norm.fit(chi[chi_mask])
        plt.hist(chi, range=(-10,10), bins=36, density=True)
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        x = np.linspace(xmin, xmax, 100)
        y = norm.pdf(x, mean, std)
        plt.vlines(mean+std, 0, ymax, color='r', linewidth=2, ls='dashed')
        plt.vlines(mean-std, 0, ymax, color='r', linewidth=2, ls='dashed')
        plt.plot(x, y, 'r-', lw=2, label='norm pdf')
        plt.title('noise check', fontsize=36)
        plt.tick_params(labelsize=20)
        plt.xlabel(r'$\chi=\frac{f_{\rm noiseless,\lambda}-f_{\rm noisy,\lambda}}{\sigma_{\lambda}}$', fontsize=36)
        plt.show()

        # sanity check - 3: vs. real object
        plt.figure(figsize=(20,6))
        fake_counts_sm, fake_ivar_sm = ivarsmooth(sobjs_fake.BOX_COUNTS, sobjs_fake.BOX_COUNTS_IVAR, 5)
        real_counts_sm, real_ivar_sm = ivarsmooth(sobjs[trace_id].BOX_COUNTS, sobjs[trace_id].BOX_COUNTS_IVAR, 5)
        plt.plot(sobjs_fake.BOX_WAVE, fake_counts_sm, drawstyle='steps-mid', color='black', label='fake')
        plt.plot(sobjs[trace_id].BOX_WAVE, real_counts_sm, drawstyle='steps-mid', color='blue', label='real')
        plt.plot(sobjs_fake.BOX_WAVE, np.sqrt(inverse(fake_ivar_sm)), drawstyle='steps-mid', color='red')
        plt.plot(sobjs_fake.BOX_WAVE, np.sqrt(inverse(real_ivar_sm)), drawstyle='steps-mid', color='orange')

        # plt.plot(sobjs_fake.BOX_WAVE, sobjs_fake.BOX_COUNTS, drawstyle='steps-mid', label='fake')
        # plt.plot(sobjs[trace_id].BOX_WAVE, sobjs[trace_id].BOX_COUNTS, drawstyle='steps-mid', color='black', label='real')
        # plt.plot(sobjs_fake.BOX_WAVE, np.sqrt(inverse(sobjs_fake.BOX_COUNTS_IVAR)), drawstyle='steps-mid', color='red')
        # plt.plot(sobjs[trace_id].BOX_WAVE, np.sqrt(inverse(sobjs[trace_id].BOX_COUNTS_IVAR)), drawstyle='steps-mid', color='orange')
        plt.ylim(np.median(sobjs_fake.BOX_COUNTS)-3*np.median(np.sqrt(inverse(sobjs_fake.BOX_COUNTS_IVAR))), 
                 np.median(sobjs_fake.BOX_COUNTS)+8*np.median(np.sqrt(inverse(sobjs_fake.BOX_COUNTS_IVAR))))
        plt.title('fake vs. real object', fontsize=40)
        plt.xlabel('Wavelength [Angstrom]', fontsize=40)
        plt.ylabel('Counts', fontsize=40)
        plt.tick_params(labelsize=20)
        plt.legend(fontsize=20)
        plt.show()
        embed()

    return img_fake, sobjs_fake, info

# hdul = io.fits_open("../resource/sensfunc/GD153_lris_long_8.7_sens.fits")
hdul = io.fits_open("../resource/sensfunc/GD153_lris_sens.fits")
sens = IRSensFunc.from_hdu(hdul)

sim_lris_path = '../resource/simulation/lris/'
targets = ['J1319+0101', 'J0901+2906', 'J1724+3718', 'J1326+0927']
trace_ids = [2, 2, 1, 2]
redshifts = [5.726, 6.093, 5.733, None]
m_Js = [21.20, 20.9, 21, 20.12]

i = 3
target = targets[i]
trace_id = trace_ids[i]
redshift = redshifts[i]
m_J = m_Js[i]

hdul = io.fits_open(os.path.join(sim_lris_path, f'{target}_tellmodel.fits'))
telluric = hdul[1].data

spec2dfile = os.path.join(sim_lris_path, f'spec2d_{target}.fits')
spec1dfile = spec2dfile.replace('spec2d', 'spec1d')

detname = DetectorContainer.get_name(det=1)
spec2DObj = spec2dobj.Spec2DObj.from_file(spec2dfile, detname, chk_version=False)
sobjs = specobjs.SpecObjs.from_fitsfile(spec1dfile, chk_version=False)
header = fits.getheader(spec1dfile)

# 1. simulate the target
# img_fake, sobjs_fake, info = simulate(sens=sens, spec2DObj=spec2DObj, sobjs=sobjs, telluric=telluric, slitloss=False, header=header, trace_id=trace_id, offset=-100, exptime=50, 
#                                       load_func=load_quasar, parse_func=parse_quasar, show_trace=False, redshift=redshift, m_J=m_J, debug=True)
img_fake, sobjs_fake, info = simulate(sens=sens, spec2DObj=spec2DObj, sobjs=sobjs, telluric=telluric, slitloss=False, header=header, trace_id=trace_id, offset=-100, exptime=300, 
                                      load_func=load_star, parse_func=parse_star, show_trace=False, m_J=m_J, debug=True)

# 2. various magnitude
# redshifts = np.arange(5.5, 7.7, 0.1)
# snr = np.zeros_like(redshifts)
# for i, redshift in enumerate(redshifts):
#     img_fake, sobjs_fake, info = simulate(sens=sens, spec2DObj=spec2DObj, sobjs=sobjs, telluric=telluric, slitloss=False, header=header, trace_id=trace_id, offset=-100, exptime=300, 
#                                           load_func=load_quasar, parse_func=parse_quasar, show_trace=False, redshift=redshift, m_J=22, debug=False, verbose=False)
#     snr[i] = info['SNR']

# 3. simulate the target
# img_fake, sobjs_fake, info = simulate(sens=sens, spec2DObj=spec2DObj, sobjs=sobjs, telluric=telluric, slitloss=False, header=header, trace_id=trace_id, offset=-100, exptime=80, 
#                                       load_func=load_quasar, parse_func=parse_quasar, show_trace=False, redshift=6.5, m_J=22, debug=True)