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

targets = ['J0901+2906', 'J1319+0101', 'J1458+1012', 'J1724+3718', 'J1111+0640', 'J1401+4542', 'J1523+2935']
redshifts = [6.1, 5.7, 5.6, 5.74, 6.0, 5.5, 5.7]
# survyes = ['UKI/VIK', 'UKI/VIK', 'PS1', 'PS1', 'UKI/VIK', 'LoTSS', 'LoTSS']
survyes = ['XDHZQSO', 'XDHZQSO', 'PS1COLOR', 'PS1COLOR', 'XDHZQSO', 'LOFAR', 'LOFAR']


spec2dfiles = ['../arxiv/LRIS_2203/LRIS_220306/reduced/all/Science/spec2d_r220306_00044-J0901+2906_OFF_LRISr_20220306T054911.654.fits',
               '../arxiv/LRIS_2203/LRIS_220305/reduced/all/Science/spec2d_r220305_00190-J1319+0101_OFF_LRISr_20220305T143022.320.fits',
               '../arxiv/LRIS_2203/LRIS_220306/reduced/J1458/Science/spec2d_r220306_00154-J1458+1012_OFF_LRISr_20220306T144525.546.fits',
               '../arxiv/LRIS_2203/LRIS_220306/reduced/all/Science/spec2d_r220306_00165-J1724+3718_OFF_LRISr_20220306T154632.448.fits',
               '../arxiv/LRIS_2204/reduced/all/Science/spec2d_r220423_00049-J1111+0640_OFF_LRISr_20220423T081759.280.fits',
               '../arxiv/LRIS_2204/reduced/all/Science/spec2d_r220423_00100-J1401+4542_OFF_LRISr_20220423T131547.837.fits',
               '../arxiv/LRIS_2204/reduced/all/Science/spec2d_r220423_00084-J1523+2935_OFF_LRISr_20220423T113251.619.fits',
               ]
spec1dfiles = ['../arxiv/LRIS_2203/LRIS_220306/reduced/all/J0901+2906/J0901+2906_coadd_tellcorr.fits',
               '../arxiv/LRIS_2203/LRIS_220305/reduced/all/J1319+0101/J1319+0101_coadd_tellcorr.fits',
               '../arxiv/LRIS_2203/LRIS_220306/reduced/J1458/J1458+1012_coadd_tellcorr.fits',
               '../arxiv/LRIS_2203/LRIS_220306/reduced/all/J1724+3718/J1724+3718_coadd_tellcorr.fits',
               '../arxiv/LRIS_2204/reduced/all/J1111+0640/J1111+0640_coadd_tellcorr.fits',
               '../arxiv/LRIS_2204/reduced/all/J1401+4542/J1401+4542_coadd_tellcorr.fits',
               '../arxiv/LRIS_2204/reduced/all/J1523+2935/J1523+2935_coadd_tellcorr.fits',
               ]

centers = [939, 931, 939, 931, 940, 939, 940]
ylims = []
for i in range(len(centers)):
    ylims.append([centers[i]-50, centers[i]+50])
xrange = (7300, 10600)

# order all the lists in the redshift order
redshifts = np.array(redshifts)
idx = np.argsort(redshifts)
targets = np.array(targets)[idx]
redshifts = redshifts[idx]
survyes = np.array(survyes)[idx]
spec2dfiles = np.array(spec2dfiles)[idx]
spec1dfiles = np.array(spec1dfiles)[idx]
ylims = np.array(ylims)[idx]

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
    if redshifts[i] > 0:
        title = targets[i] + r', $z$ = ' + str(redshifts[i]) + f', {survyes[i]}'
    else:
        title = targets[i] + f', {survyes[i]}'
    ax, wl, tell = plot_spec1d(title, spec1dfiles[i], idx=1, axis=ax2, smooth_window=5, 
                               template=False, telluric=True, qso=True, plot_telluric=False)
    # interpolate telluric, and apply on a grid
    from scipy.interpolate import interp1d
    if i == 0:
        wl_grid = np.arange(wl[0], wl[-1] , 0.1)
    tell_interp = interp1d(wl, tell, fill_value="extrapolate")
    tell_grid = tell_interp(wl_grid)
    tell_grid[tell_grid < 0] = 0
    tell_grid[tell_grid > 1] = 1
    telluric_all.append(tell_grid)

    if targets[i] == 'J1523+2935' or targets[i] == 'J1724+3718':
        ax2.set_ylim(ax2.get_ylim()[0], 1.2 * ax2.get_ylim()[1])
    ax2.legend(loc="upper right", frameon=True, fontsize=12)
    ax2.tick_params(axis='both', which='major', labelsize=20, width=1, size=6)
    ax2.tick_params(axis='both', which='minor', labelsize=20, width=1, size=3)
    ax2.set_xlim(xrange[0], xrange[1])

    fig.subplots_adjust(hspace=0)

axs[-1].plot(wl_grid, np.mean(np.array(telluric_all), axis=0), color="grey", zorder=1, alpha=0.8)
axs[-1].set_ylim(0., 1.5)
axs[-1].set_xlim(xrange[0], xrange[1])
axs[-1].set_xlabel(r"Wavelength ($\mathring{A}$)", fontsize=20)
axs[-1].tick_params(axis='both', which='major', labelsize=20, width=1, size=6)
axs[-1].tick_params(axis='both', which='minor', labelsize=20, width=1, size=3)
# text on lower right
axs[-1].text(0.975, 0.1, "average telluric model", horizontalalignment='right', verticalalignment='bottom', transform=axs[-1].transAxes, fontsize=15)

axs[len(targets)].set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\mathring{A}^{-1}$)", fontsize=20)

if args.show:
    plt.show()
else:
    plt.savefig(f'qso_1d2d.pdf')