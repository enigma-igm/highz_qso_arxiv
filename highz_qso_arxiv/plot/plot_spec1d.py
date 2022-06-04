from highz_qso_arxiv.util import ivarsmooth, inverse

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
plt.style.use('science')

import warnings
warnings.filterwarnings("ignore")

from IPython import embed

ARXIV_PATH = "../arxiv/"

def plot_spec1d(name, fits_file, idx, axis, smooth_window=5, telluric=False, telluric_fits_file=""):
    """Plot single spectrum to axis

    Args:
        name (str): _description_
        fits_file (str): _description_
        idx (int): _description_
        axis (matplotlib.axes): _description_
        smooth_window (int, optional): _description_. Defaults to 5.
        telluric (bool, optional): _description_. Defaults to False.
        telluric_fits_file (str, optional): _description_. Defaults to "".

    Returns:
        _type_: _description_
    """
    hdul = fits.open(fits_file)
    data = Table(hdul[idx].data)
    try:
        wave, flux, flux_ivar = data["wave"], data["flux"], data["ivar"]
    except KeyError:
        wave, flux, flux_ivar = data["OPT_WAVE"], data["OPT_FLAM"], data["OPT_FLAM_IVAR"]
    
    # TODO: mask
    mask = wave>0
    flux_sm, flux_ivar_sm = ivarsmooth(flux, flux_ivar, smooth_window)
    axis.plot(wave[mask], flux_sm[mask], label=name, color="black", lw=1.5)
    axis.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm[mask])), 
              lw=1, color="red", alpha=0.6)
    axis.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
    axis.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
    axis.set_xlim(np.min(wave), np.max(wave))

    # always want to include the noise vector
    ymin = np.mean(flux_sm)-2*np.std(flux_sm)
    if ymin > 0: ymin = 0
    axis.set_ylim(ymin,np.mean(flux_sm)+2*np.std(flux_sm))
    axis.legend(loc="upper right", frameon=True)
    if telluric:
        hdul_tell = fits.open(telluric_fits_file)
        data = Table(hdul_tell[1].data)
        wave, flux = data["WAVE"].value[0], data["TELLURIC"].value[0]
        ax2 = axis.twinx()
        ax2.plot(wave[mask], flux[mask], color="navy", zorder=1, alpha=0.4)
        ax2.set_ylim(-2,1.)
        ax2.set_yticks([0,0.5,1])

    return axis

def plot_single(name, fits_file, idx, smooth_window=5, telluric=False, telluric_fits_file="", plot=True, save_file=""):
    """Plot single spectrum given fits file and other parameters

    Args:
        name (_type_): _description_
        fits_file (_type_): _description_
        idx (_type_): _description_
        smooth_window (int, optional): _description_. Defaults to 5.
        telluric (bool, optional): _description_. Defaults to False.
        telluric_fits_file (str, optional): _description_. Defaults to "".
        plot (bool, optional): _description_. Defaults to True.
        save_file (str, optional): _description_. Defaults to "".

    Returns:
        _type_: _description_
    """
    fig, ax = plt.subplots(figsize=(20,6))
    plot_spec1d(name, fits_file, idx, ax, smooth_window, telluric, telluric_fits_file)
    if plot:
        plt.show()
    if save_file:
        fig.savefig(save_file)
    return fig, ax

def plot_series(name_list, fits_list, idx_list, smooth_window=5, plot=True, save_file=""):
    """Plot a series of spectrum given fits files and other parameters

    Args:
        name_list (_type_): _description_
        fits_list (_type_): _description_
        idx_list (_type_): _description_
        smooth_window (int, optional): _description_. Defaults to 5.
        plot (bool, optional): _description_. Defaults to True.
        save_file (str, optional): _description_. Defaults to "".

    Returns:
        _type_: _description_
    """
    num = len(fits_list)
    fig, axs = plt.subplots(num, 1, figsize=(12,3*num))
    for i, ax in enumerate(axs):
        # TODO: telluric for each plot
        plot_spec1d(name_list[i], fits_list[i], idx_list[i], ax, smooth_window)
    fig.tight_layout()
    if plot:
        plt.show()
    if save_file:
        fig.savefig(save_file)
    return fig, axs

if __name__ == "__main__":
    name = "J2212"
    fits_file = ARXIV_PATH + "J2212/J2212+2040_LRIS_coadd1d_tellcorr.fits"
    plot_single(name, fits_file, 1)