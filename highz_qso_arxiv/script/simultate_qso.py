import matplotlib.pyplot as plt
from pypeit.core.flux_calib import Flam_to_Nlam
from pypeit.core.wavecal import wvutils
from astropy.io import ascii, fits
from scipy import interpolate

from IPython import embed

def simulate_qso_Npix(redshift, exptime):
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
    wl_lya = 1215.67 * (1 + redshift)
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

    fig, ax = plt.subplots(figsize=(12,6))
    ax.plot(wl_obs[mask], Npix)
    ax.set_xlim(wl_sens[wl_gpm][0], wl_sens[wl_gpm][-1])
    ax.set_xlabel(r"Wavelength [$\AA$]")
    ax.set_ylabel(r"$\rm N_{pix}$ [$\mathrm{photons/pixel}$]")
    plt.show()
    
    return

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

if __name__ == "__main__":
    simulate_qso_Npix(6.5, 20)
    # plot_qso(6) 