from ..util import ivarsmooth, inverse

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits, ascii

import warnings
warnings.filterwarnings("ignore")

from IPython import embed

ARXIV_PATH = "/Volumes/Extreme SSD/highz_qso_arxiv/highz_qso_arxiv/arxiv/"
RESOURCE_PATH = "/Volumes/Extreme SSD/highz_qso_arxiv/highz_qso_arxiv/resource/"

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
    axis.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
    axis.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
    axis.set_xlim(np.min(wave[mask]), np.max(wave[mask]))

    # always want to include the noise vector
    
    # TODO: dynamically determine the ylim

    mask = flux_std < .5
    # flux_oversm, _ = ivarsmooth(flux, flux_ivar, 11)
    ymin = np.mean(flux_sm[mask])-2*np.std(flux_sm[mask])
    ymax = np.mean(flux_sm[mask])+2*np.std(flux_sm[mask])
    if ymin > 0: ymin = 0
    # if ymax < np.max(flux_sm[mask]): ymax =  np.mean(flux_sm)+3*np.std(flux_sm)
    axis.set_ylim(ymin, ymax)

    if template:
        from scipy.interpolate import interp1d
        templates = ["L0", "L0.5", "L1", "L1.5", "L2", "L3", "L5", "L8", "L6",
                     "M4.5", "M5", "M6", "M7", "M8", "M9", "M9.5",
                     "T0"]
        chi_sq = np.inf
        for template in templates:
            # TODO: bin multiple templates for each type
            tab = ascii.read(RESOURCE_PATH+f"dwarf/keck_lris_{template}_1.dat")
            wave_template, flux_template = tab["col1"], tab["col2"]/1e-17
            template_interp_func = interp1d(wave_template, flux_template, kind="cubic")
            mask_for_scale = wave < min(np.max(wave_template), np.max(wave))

            # interpolate the template to the same wavelength as the data
            flux_template_interp = template_interp_func(wave[mask_for_scale])

            # calculate the scale factor
            # d(chi_sq)/d(scale_factor) = 0
            scale = np.sum(flux_sm[mask_for_scale]*flux_template_interp/flux_std[mask_for_scale]**2) / \
                    np.sum(flux_template_interp**2/flux_std[mask_for_scale]**2)
            chi_sq_new = np.sum((flux_sm[mask_for_scale]-scale*flux_template_interp)**2/flux_std[mask_for_scale]**2)

            # choose the template with the lowest chi_sq
            if chi_sq_new < chi_sq:
                chi_sq = chi_sq_new
                best_template = (wave[mask_for_scale], scale*flux_template_interp, template, chi_sq)
        axis.plot(best_template[0], best_template[1], alpha=0.8, label=best_template[2])  
                #   label=best_template[2]+"\n"+r"$\chi^2=$"+str(round(best_template[3],2)))

    if telluric:
        flux =data["telluric"]
        ax2 = axis.twinx()
        ax2.plot(wave[mask], flux[mask], color="navy", zorder=1, alpha=0.4)
        ax2.set_ylim(-2,1.)
        ax2.set_yticks([0,0.5,1])

    axis.legend(loc="upper right", frameon=True)
    return axis

def plot_single(name, fits_file, idx, smooth_window=5, telluric=False, display=True, save_file=""):
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
    plot_spec1d(name, fits_file, idx, ax, smooth_window, telluric)
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
    fig, axs = plt.subplots(num, 1, figsize=(12,3*num))
    for i, ax in enumerate(axs):
        if template_list is None: template = True
        else: template = template_list[i]
        if telluric_list is None: telluric = False
        else: telluric = telluric_list[i]
        plot_spec1d(name_list[i], fits_list[i], idx_list[i], ax, smooth_window, template=template, telluric=telluric)
    fig.tight_layout()
    if display:
        plt.show()
    if save_file:
        fig.savefig(save_file)
    return fig, axs

if __name__ == "__main__":
    name = "J2212"
    fits_file = ARXIV_PATH + "J2212/J2212+2040_LRIS_coadd1d_tellcorr.fits"
    plot_single(name, fits_file, 1)