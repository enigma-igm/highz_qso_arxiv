from re import template
from highz_qso_arxiv.plot import plot_spec1d
from highz_qso_arxiv.util import ivarsmooth, inverse
from highz_qso_arxiv.util.spec1dutil import rescale
from matplotlib.ticker import MaxNLocator

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits, ascii

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
        templates = ["L0.5", "L6",   "M7",
                     "L1",   "L7.5", "M8",
                     "L2",   "L8",   "M9.5",
                     "L3.5", "M4.5", "M9",
                     "L3",   "M5",   "T2",
                     "L4.5", "M6.5", "T4.5",
                     "L5",   "M6"]

        chi_sq = np.inf
        for template in templates:
            # tab = ascii.read(os.path.join(RESOURCE_PATH, f"dwarf/keck_lris_{template}_1.dat"))
            try:
                tab = ascii.read(os.path.join(RESOURCE_PATH, f"dwarf/irtf_{template}_1.dat"))
            except:
                print(template)
            # wave_template, flux_template = tab["col1"], tab["col2"]/1e-17
            wave_template, flux_template = tab["col1"]*1e4, tab["col2"]

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

    ax.set_xlabel(r"wavelength ($\AA$)", fontsize=20)
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

# # LRIS-2201
# name_list = ['J0759+2811', 'J0739+2328', 'J0820+2412', 'J0831+2558',
#              'J0849+2601', 'J0936+3244', 'J0936+3346', 'J0938+3332',
#              'J1030-0031', 'J1003-2610', 'J1145+0933', 'J1146+0039',
#              'J1206-0051', 'J1241-0134', 'J1246-3045', 'J1256-0306',
#              'J1305-1549', 'J1326+0927', 'J1209+0135']
# fits_list = [f"../arxiv/LRIS_2201/reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in name_list]
# fits_list[10] = f"../arxiv/LRIS_2201/reduced/all/coadd2d/{name_list[10]}_coadd_tellcorr.fits"
# fits_list[12] = f"../arxiv/LRIS_2201/reduced/all/coadd2d/{name_list[12]}_coadd_tellcorr.fits"
# fits_list[15] = f"../arxiv/LRIS_2201/reduced/all/coadd2d/{name_list[15]}_coadd_tellcorr.fits"
# idx_list = [1 for i in range(len(name_list))]

# # targets_0305 = ["J1223+0114","J1250+0347", "J1319+0101"]
# # fits_list = [f"../arxiv/LRIS_2203/LRIS_220305/reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in targets_0305]
# # idx_list = [1 for i in range(len(targets_0305))]

# # a list of True with the same length as fits_list
# template_list = [True for i in range(len(fits_list))]
# template_list[2] = False
# telluric_list = [True for i in range(len(fits_list))]
# plot_series(name_list, fits_list, idx_list, smooth_window=3, 
#             template_list=template_list, telluric_list=telluric_list, display=False, save_file="LRIS_2201.pdf")

# # LRIS-2203
# targets_0305 = ["J1100+0203", "J1143-0248", "J1200+0112", "J1223+0114",
#                 "J1250+0347", "J1319+0101", "J1327-0207", "J1335+0103",
#                 "J1355-0044", "J1355-0111", "J1356+0050", "J1358-0118",
#                 "J1422+0307"]
# fits_list_0305 = [f"../arxiv/LRIS_2203/LRIS_220305/reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in targets_0305]
# fits_list_0305[0] = "../arxiv/LRIS_2203/LRIS_220305/reduced/all/coadd2d/J1100+0203_coadd_tellcorr.fits"
# fits_list_0305[2] = "../arxiv/LRIS_2203/LRIS_220305/reduced/all/coadd2d/J1200+0112_coadd_tellcorr.fits"

# # a list of True with the same length as fits_list
# template_list_0305 = [True for i in range(len(fits_list_0305))]
# template_list_0305[5] = False
# telluric_list_0305 = [True for i in range(len(fits_list_0305))]
# plot_series(targets_0305, fits_list_0305, [1 for i in range(len(fits_list_0305))], smooth_window=3,
#             template_list=template_list_0305, telluric_list=telluric_list_0305, display=False, save_file="LRIS_2203_0305.pdf") 

# targets_0306 = ["J0810+2352", "J0850+0146", "J0901+2906", "J0911+0022",
#                 "J0947+0111", "J1153-2239", "J1233-2807", "J1412-2757",
#                 "J1425+0004", "J1427+0202", "J1458+1012", "J1510-0144", 
#                 "J1514-0121", "J1522+0051", "J1534+0046", "J1655-0051", "J1724+3718"]
# fits_list_0306 = [f"../arxiv/LRIS_2203/LRIS_220306/reduced/all/{nm}/{nm}_coadd_tellcorr.fits" for nm in targets_0306]
# fits_list_0306[9] = "../arxiv/LRIS_2203/LRIS_220306/reduced/all/coadd2d/J1427+0202_coadd_tellcorr.fits"

# template_list_0306 = [True for i in range(len(fits_list_0306))]
# template_list_0306[2] = False
# template_list_0306[10] = False
# template_list_0306[-1] = False
# telluric_list_0306 = [True for i in range(len(fits_list_0306))]
# plot_series(targets_0306, fits_list_0306, [1 for i in range(len(fits_list_0306))], smooth_window=3,
#             template_list=template_list_0306, telluric_list=telluric_list_0306, display=False, save_file="LRIS_2203_0306.pdf")

# targets = targets_0305 + targets_0306
# fits_list = fits_list_0305 + fits_list_0306
# template_list = template_list_0305 + template_list_0306
# telluric_list = telluric_list_0305 + telluric_list_0306
# idx_list = [1 for i in range(len(targets))]
# plot_series(targets, fits_list, idx_list, smooth_window=3,
#             template_list=template_list, telluric_list=telluric_list, display=False, save_file="LRIS_2203.pdf")

# LRIS-2204
name_list = ["J0739+2645", "J0954-0117", "J1111+0640",
             "J1238+0234", "J1253-0200", "J1317-0302",
             "J1340+1422", "J1349+0805", "J1351+0128", "J1352+0002",
             "J1438-0103", "J1426+0148", "J1429+0219",
             "J1434+0240", "J1523+2935", "J1625+2414", "J1635+5940",
             "J1441+0149", "J1437-0151", "J1450+0228", "J1401+4542",
             "J1535+6146", "J1508-0109", "J1513-0059", "J1532-0155",
             "J2212+2040"]
fits_list = [f"../arxiv/LRIS_2204/reduced/all/{nm}/{nm}_coadd.fits" for nm in name_list]
fits_list[-1] = "../arxiv/LRIS_2204/reduced/J2212/J2212+2040_coadd.fits"
idx_list = [1 for i in range(len(name_list))]
template_list = [True for i in range(len(fits_list))]
telluric_list = [False for i in range(len(fits_list))]
plot_series(name_list, fits_list, idx_list, template_list=template_list, telluric_list=telluric_list, 
            display=False, save_file="LRIS_2204.pdf")