import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval, ImageNormalize, SqrtStretch

from highz_qso_arxiv.plot import plot_spec1d

from pypeit import spec2dobj
from pypeit.images.detector_container import DetectorContainer

from IPython import embed

parser = argparse.ArgumentParser()
parser.add_argument('--show', action='store_true')
args = parser.parse_args()

targets = ['J0739+2328', 'J1209+0135', 'J1238+0234', 'J1327-0207']
types = ['M9.5', 'L0.5', 'UNQ', 'UNQ']
survyes = ['XDHZQSO', 'XDHZQSO', 'XDHZQSO', 'XDHZQSO']

spec2dfiles = ['../arxiv/LRIS_2201/reduced/all/Science/spec2d_r220127_00086-J0739+2328_OFF_LRISr_20220127T102658.646.fits',
               '../arxiv/LRIS_2201/reduced/all/Science/spec2d_r220127_00123-J1209+0135_OFF_LRISr_20220127T145655.277.fits',
               '../arxiv/LRIS_2204/reduced/all/Science/spec2d_r220423_00053-J1238+0234_OFF_LRISr_20220423T084038.957.fits',
               '../arxiv/LRIS_2203/LRIS_220305/reduced/all/Science/spec2d_r220305_00192-J1327-0207_OFF_LRISr_20220305T144445.197.fits']
spec1dfiles = ['../arxiv/LRIS_2201/reduced/all//J0739+2328/J0739+2328_coadd_tellcorr.fits',
               '../arxiv/LRIS_2201/reduced/all//J1209+0135/J1209+0135_coadd_tellcorr.fits',
               '../arxiv/LRIS_2204/reduced/all//J1238+0234/J1238+0234_coadd_tellcorr.fits',
               '../arxiv/LRIS_2203/LRIS_220305/reduced/all//J1327-0207/J1327-0207_coadd_tellcorr.fits']

centers = [941, 941, 938, 943]
ylims = []
for i in range(len(centers)):
    ylims.append([centers[i]-50, centers[i]+50])
xrange = (7300, 10600)

height_ratios = []
for _ in range(len(targets)):
    height_ratios += [1, 2]
height_ratios += [1]

fig, axs = plt.subplots(2*len(targets)+1, 1, figsize=(15,2*len(targets)), gridspec_kw={'height_ratios': height_ratios})

telluric_all = []
for i in range(len(targets)):
    detname = DetectorContainer.get_name(det=1)
    spec2DObj = spec2dobj.Spec2DObj.from_file(spec2dfiles[i], detname, chk_version=False)

    gpm = (spec2DObj.bpmmask == 0)
    # image = (spec2DObj.sciimg - spec2DObj.skymodel) * np.sqrt(spec2DObj.ivarmodel) * gpm
    image = (spec2DObj.sciimg - spec2DObj.skymodel) * gpm

    from astropy.stats import sigma_clipped_stats
    mean, med, sigma = sigma_clipped_stats(image[spec2DObj.bpmmask == 0], sigma_lower=5.0,
                                        sigma_upper=5.0)
    cut_min = mean - 1.0 * sigma
    cut_max = mean + 4.0 * sigma

    norm = ImageNormalize(image, interval=ZScaleInterval(), stretch=SqrtStretch(), vmin=cut_min, vmax=cut_max)

    ax1 = axs[i*2]
    ax2 = axs[i*2+1]

    ax1.imshow(image.T, origin='lower', norm=norm, cmap='gray', aspect = "auto")
    wv = spec2DObj.waveimg[:,centers[i]]
    idx = np.where((wv > xrange[0]) & (wv < xrange[1]))[0]
    ax1.set_xlim(idx[0], idx[-1])
    ax1.set_ylim(ylims[i])
    ax1.set_xticks([])
    ax1.set_yticks([])
    title = f'{targets[i]}, {types[i]}, {survyes[i]}'
    if types[i] == 'UNQ':
        template = False
    else:
        template = True
    ax, wl, tell, chi_sq = plot_spec1d(title, spec1dfiles[i], idx=1, axis=ax2, smooth_window=5, 
                                       template=True, telluric=True, qso=False, plot_telluric=False)

    # interpolate telluric, and apply on a grid
    from scipy.interpolate import interp1d
    if i == 0:
        wl_grid = np.arange(wl[0], wl[-1] , 0.1)
    tell_interp = interp1d(wl, tell, fill_value="extrapolate")
    tell_grid = tell_interp(wl_grid)
    tell_grid[tell_grid < 0] = 0
    tell_grid[tell_grid > 1] = 1
    telluric_all.append(tell_grid)
 
    if targets[i] == 'J0739+2328' or targets[i] == 'J1209+0135':
        ax2.set_ylim(ax2.get_ylim()[0], 1.5 * ax2.get_ylim()[1])
    elif targets[i] == 'J1238+0234':
        ax2.set_ylim(ax2.get_ylim()[0], 1.2 * ax2.get_ylim()[1])
    ax2.legend(loc="upper left", frameon=True, fontsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=20, width=1, size=6)
    ax2.tick_params(axis='both', which='minor', labelsize=20, width=1, size=3)
    ax2.tick_params(axis='both', which='major', labelsize=18)
    ax2.set_xlim(xrange[0], xrange[1])
    
    fig.subplots_adjust(hspace=0)

axs[-1].plot(wl_grid, np.mean(np.array(telluric_all), axis=0), color="grey", zorder=1, alpha=0.8)
axs[-1].set_ylim(0., 1.5)
axs[-1].set_xlim(xrange[0], xrange[1])
axs[-1].set_xlabel(r"Wavelength ($\mathring{A}$)", fontsize=20)
axs[-1].tick_params(axis='both', which='major', labelsize=20, width=1, size=6)
axs[-1].tick_params(axis='both', which='minor', labelsize=20, width=1, size=3)
axs[-1].text(0.975, 0.1, "average telluric model", horizontalalignment='right', verticalalignment='bottom', transform=axs[-1].transAxes, fontsize=15)

axs[len(targets)].set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\mathring{A}^{-1}$)", fontsize=20)
axs[len(targets)].yaxis.set_label_coords(-0.04, 0)

if args.show:
    plt.show()
else:
    plt.savefig(f'noqso_1d2d.pdf')