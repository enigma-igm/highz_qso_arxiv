from ..util.spec1dutil import rescale
from ..util import ivarsmooth, inverse, get_project_root

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits, ascii
from matplotlib.ticker import MaxNLocator

import warnings
warnings.filterwarnings("ignore")

from IPython import embed

path = get_project_root()
ARXIV_PATH = os.path.join(path, "arxiv")
RESOURCE_PATH = os.path.join(path, "resource")

def plot_template_dat(model='qso', redshift=7., star_type='L0', display=True):
    """Available template models:
        qso: QSO template (Selsing+2015)
        star: Dwarf template (L0, L0.5, L1, L1.5, L2, L3, L5, L6, L8, 
                              M4.5, M5, M6, M7, M8, M9, M9.5, T0)
    Args:
        model (str, optional): _description_. Defaults to 'qso'.
        redshift (_type_, optional): _description_. Defaults to 7..
        star_type (str, optional): _description_. Defaults to 'L0'.
        display (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """
    if model == 'qso':
        dat = ascii.read(os.path.join(RESOURCE_PATH, "Selsing2015.dat"))
        wl_rest = dat["col1"]
        wl_obs = wl_rest * (1 + redshift)
        flux = dat["col2"] # in 1e-17 erg/s/cm2/A
        flux_err = dat["col3"]
        wl_lya = 1215.67 * (1 + redshift)
        trough = wl_obs < wl_lya
        flux[trough] = 0
        label = rf"$z={redshift}$; Selsing 2015"
    elif model == 'star':
        dat = ascii.read(os.path.join(RESOURCE_PATH, f"dwarf/keck_lris_{star_type}_1.dat"))
        wl_obs = dat["col1"]
        flux = dat["col2"] * 1e17 # in 1e-17 erg/s/cm2/A
        label = star_type
    fig, ax = plt.subplots(figsize=(12,6))
    ax.plot(wl_obs, flux, label=label, color="black")
    # ax.fill_between(wl_obs, flux-flux_err, flux+flux_err, color="black", alpha=0.2)
    ax.set_xlim(7300, 10500)
    ax.legend(loc="upper left")
    ax.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
    ax.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)

    # ymin, ymax = ax.get_ylim()
    # ax.vlines(wl_lya, ymin, ymax)
    if display:
        plt.show()
    return fig, ax

def plot_spec1d(name, fits_file, idx, axis, smooth_window=5, template=True, telluric=False):
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
    mask = wave>1.
    flux_sm, flux_ivar_sm = ivarsmooth(flux, flux_ivar, smooth_window)
    flux_std = inverse(np.sqrt(flux_ivar_sm[mask]))
    axis.plot(wave[mask], flux_sm[mask], label=name, color="black", lw=1.5)
    axis.plot(wave[mask], flux_std, lw=1, color="red", alpha=0.6)
    xmin = np.min(wave[mask])
    xmax = np.max(wave[mask])
    axis.hlines(0, xmin, xmax, color="cyan", alpha=0.8, lw=1, ls='dashed')
    axis.set_xlim(xmin, xmax)

    # always want to include the noise vector
    
    # TODO: dynamically determine the ylim

    mask = flux_std < .5
    # flux_oversm, _ = ivarsmooth(flux, flux_ivar, 11)
    ymin = np.mean(flux_sm[mask])-2*np.std(flux_sm[mask])
    ymax = np.mean(flux_sm[mask])+2*np.std(flux_sm[mask])
    if ymin > 0: ymin = 0
    # if ymax < np.max(flux_sm[mask]): ymax =  np.mean(flux_sm)+3*np.std(flux_sm)
    axis.set_ylim(ymin, ymax)
    axis.yaxis.set_major_locator(MaxNLocator(2, prune='upper'))
    yt = axis.get_yticks()
    if yt[-1] + 0.1 > ymax:
        axis.set_ylim(ymin, ymax+0.2)
        axis.set_yticks(yt[1:])
    # axis.set_yticks([])

    if template:
        templates = ["L0.5", "L5",   "M7",
                     "L1",   "M4.5", "M8",
                     "L2",   "M5",   "M9.5",
                     "L3",   "M6",   "M9", "L6", "L8", "T"]
        chi_sq = np.inf
        for template in templates:
            try:
                tab = ascii.read(os.path.join(RESOURCE_PATH, f"dwarf/combine_lris_{template}.dat"))
            except:
                print(template)
            wave_template, flux_template = tab["wave"], tab["flux"]
            flux_template_scaled, chi_sq_new = rescale(wave_template, flux_template, 
                                                       wave, flux_sm, flux_std)

            # choose the template with the lowest chi_sq
            # print(template, chi_sq_new)
            if chi_sq_new < chi_sq:
                chi_sq = chi_sq_new
                best_template = (wave_template, flux_template_scaled, template, chi_sq)
        axis.plot(best_template[0], best_template[1], alpha=0.8, label=best_template[2])  
                #   label=best_template[2]+"\n"+r"$\chi^2=$"+str(round(best_template[3],2)))

    if telluric:
        flux =data["telluric"]
        ax2 = axis.twinx()
        ax2.plot(wave[mask], flux[mask], color="navy", zorder=1, alpha=0.4)
        ax2.set_ylim(-2,1.)
        ax2.set_yticks([])
    axis.tick_params(axis='both', which='major', labelsize=10)
    axis.legend(loc="upper right", frameon=True, fontsize=6)
    return axis

def plot_single(name, fits_file, idx, smooth_window=5, template=True, telluric=False, display=True, save_file=""):
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
    plot_spec1d(name, fits_file, idx, ax, smooth_window, template=template, telluric=telluric)
    if display:
        plt.show()
    if save_file:
        fig.savefig(save_file)
    return fig, ax

def plot_series(name_list, fits_list, idx_list, smooth_window=5, template_list=None, telluric_list=None, display=True, save_file=""):
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
    assert len(name_list) == len(fits_list) == len(idx_list)
    if template_list is not None:
        assert len(template_list) == len(name_list)
    if telluric_list is not None:
        assert len(telluric_list) == len(name_list)

    num = len(fits_list)
    fig, axs = plt.subplots(num, sharex=True, figsize=(10,1*num))
    for i, ax in enumerate(axs):
        if template_list is None: template = True
        else: template = template_list[i]
        if telluric_list is None: telluric = False
        else: telluric = telluric_list[i]
        ax = plot_spec1d(name_list[i], fits_list[i], idx_list[i], ax, smooth_window, template=template, telluric=telluric)
    ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=20)
    mid_idx = int(len(name_list)/2)

    axs[mid_idx].set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=20)
    axs[mid_idx].yaxis.set_label_coords(-0.06, 0.5)

    # fig.tight_layout()
    fig.subplots_adjust(hspace=0)

    if display:
        plt.show()
    if save_file:
        fig.savefig(save_file)
    return fig, axs

if __name__ == "__main__":
    name = "J2212"
    fits_file = ARXIV_PATH + "J2212/J2212+2040_LRIS_coadd1d_tellcorr.fits"
    plot_single(name, fits_file, 1)