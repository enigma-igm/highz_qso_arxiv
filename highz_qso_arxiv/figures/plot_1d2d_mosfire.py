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

targets = ['J0426+0221', 'J0220+3458', 'J2346+2545', 'J0251+1724']
types = ['INCONCLUSIVE', 'INCONCLUSIVE', 'UNQ', 'UNQ']
survyes = ['Z7DROPOUT', 'LOFAR', 'Z7DROPOUT', 'Z7DROPOUT']

spec2dfiles = ['../arxiv/MOSFIRE_1911/reduced/all/coadd2d/Science_coadd/spec2d_MF.20191118.45957-MF.20191118.46603-J0426+0221.fits',
               '../arxiv/MOSFIRE_2010/reduced/all/coadd2d/Science_coadd/spec2d_MF.20201022.50645-MF.20201022.51567-lo.fits', 
               '../arxiv/MOSFIRE_2010/reduced/all/coadd2d/Science_coadd/spec2d_MF.20201022.33150-MF.20201022.34069-z7.fits',
               '../arxiv/MOSFIRE_2010/reduced/all/coadd2d/Science_coadd/spec2d_MF.20201022.46190-MF.20201022.47109-z7.fits']
spec1dfiles = ['../arxiv/MOSFIRE_1911/reduced/all/coadd2d/J0426+0221_coadd_tellcorr.fits',
               '../arxiv/MOSFIRE_2010/reduced/all/coadd2d/J0220+3458_coadd_tellcorr.fits',
               '../arxiv/MOSFIRE_2010/reduced/all/coadd2d/J2346+2545_coadd_tellcorr.fits',
               '../arxiv/MOSFIRE_2010/reduced/all/coadd2d/J0251+1724_coadd_tellcorr.fits']

centers = [1064, 1058, 1058, 1059]
ylims = []
for i in range(len(centers)):
    ylims.append([centers[i]-50, centers[i]+50])
xrange = (9650, 11200)

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
    ax2, wl, tell = plot_spec1d(title, spec1dfiles[i], idx=1, axis=ax2, smooth_window=5, 
                                template=False, telluric=True, qso=False, plot_telluric=False)

    # interpolate telluric, and apply on a grid
    from scipy.interpolate import interp1d
    if i == 0:
        wl_grid = np.arange(wl[0], wl[-1] , 0.1)
    tell_interp = interp1d(wl, tell, fill_value="extrapolate")
    tell_grid = tell_interp(wl_grid)
    tell_grid[tell_grid < 0] = 0
    tell_grid[tell_grid > 1] = 1
    telluric_all.append(tell_grid)
 
    if targets[i] == 'J0426+0221':
        ax2.set_ylim(-0.1, 0.45)
    elif targets[i] == 'J0220+3458':
        ax2.set_ylim(-0.1, 0.6)
    elif targets[i] == 'J2346+2545':
        ax2.set_ylim(-0.1, 0.6)

    ax2.legend(loc="upper right", frameon=True, fontsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=20, width=1, size=6)
    ax2.tick_params(axis='both', which='minor', labelsize=20, width=1, size=3)
    ax2.tick_params(axis='both', which='major', labelsize=18)
    ax2.set_xticks([])
    ax2.set_xlim(xrange[0], xrange[1])

    fig.subplots_adjust(hspace=0)

axs[-1].plot(wl_grid, np.mean(np.array(telluric_all), axis=0), color="grey", zorder=1, alpha=0.8)
axs[-1].set_ylim(0, 1.5)
axs[-1].set_xlim(xrange[0], xrange[1])
axs[-1].set_xlabel(r"Wavelength ($\mathring{A}$)", fontsize=20)
axs[-1].tick_params(axis='both', which='major', labelsize=20, width=1, size=6)
axs[-1].tick_params(axis='both', which='minor', labelsize=20, width=1, size=3)
axs[-1].text(0.975, 0.1, "average telluric model", horizontalalignment='right', verticalalignment='bottom', transform=axs[-1].transAxes, fontsize=15)

axs[len(targets)].set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\mathring{A}^{-1}$)", fontsize=20)
# move the ylabel down a bit
axs[len(targets)].yaxis.set_label_coords(-0.04, 0)

if args.show:
    plt.show()
else:
    plt.savefig(f'mosfire_1d2d.pdf')