import numpy as np
from highz_qso_arxiv.plot import plot_series, plot_spec1d
import matplotlib.pyplot as plt

name_all = []
fits_all = []

# NIRES-1903

name_list = ["J0953-0853", "J1638+5412", "J0756+5744", "J0800+3034"]
types = ["STAR", "STAR", "UNQ", "UNQ"]
survyes = ["Z7DROPOUT", "Z7DROPOUT", "Z7DROPOUT", "Z7DROPOUT"]

fits_list = [f"../arxiv/NIRES_1903/reduced/all/coadd2d/{nm}_coadd.fits" for nm in name_list]

name_all += name_list
fits_all += fits_list

idx_all = [1, 1, 1, 1]
template_all = [True, True, False, False]
qso_all = [False, False, False, False]

fits_all = [f[:-5] + '_tellcorr' + f[-5:] for f in fits_all]
fig, axs = plt.subplots(len(name_list)+1, 1, figsize=(15,1.5*len(name_list)), 
                        gridspec_kw={'height_ratios': [2,2,2,2,1]})

xrange = (9600, 24400)

telluric_all = []
for i in range(len(name_list)):
    ax = axs[i]
    template = template_all[i]
    title = f'{name_list[i]}, {types[i]}, {survyes[i]}'

    if template:
        ax, wl, tell, chi_sq = plot_spec1d(title, fits_all[i], 1, ax, smooth_window=9,
                                           template=template, telluric=True, qso=False, 
                                           plot_telluric=False)
    else:
        ax, wl, tell = plot_spec1d(title, fits_all[i], 1, ax, smooth_window=9,
                                   template=template, telluric=True, qso=False, 
                                   plot_telluric=False)
    # interpolate telluric, and apply on a grid
    from scipy.interpolate import interp1d
    if i == 0:
        wl_grid = np.arange(wl[0], wl[-1] , 0.1)
    tell_interp = interp1d(wl, tell, fill_value="extrapolate")
    tell_grid = tell_interp(wl_grid)
    tell_grid[tell_grid < 0] = 0
    tell_grid[tell_grid > 1] = 1
    telluric_all.append(tell_grid)

    ax.axvspan(13500, 14250, 0, 0.99, alpha=0.95, color='white', zorder=4)
    ax.axvspan(18000, 19500, 0, 0.99, alpha=0.95, color='white', zorder=4)

    ax.set_xlim(xrange[0], xrange[1])
    ax.legend(loc="upper right", frameon=True, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=16, width=1, size=4)
    ax.tick_params(axis='both', which='minor', labelsize=16, width=1, size=2)

    ax.set_ylim(-0.1, ax.get_ylim()[1] * 1.3)
    if name_list[i] == "J0756+5744":
        ax.set_ylim(-0.1, ax.get_ylim()[1] * 1.5)
        # ax.set_yticks([0, 1.0])
    elif name_list[i] == "J0800+3034":
        ax.set_ylim(-0.5, ax.get_ylim()[1] * 1.5)

axs[-1].set_ylim(0, 1.5)
axs[-1].set_yticks([0, 1])
axs[-1].set_xlim(xrange[0], xrange[1])
axs[-1].plot(wl_grid, np.mean(np.array(telluric_all), axis=0), color="grey", zorder=1, alpha=0.8)
axs[-1].set_xlabel(r"Wavelength ($\mathring{A}$)", fontsize=18)
axs[-1].tick_params(axis='both', which='major', labelsize=16, width=1, size=4)
axs[-1].tick_params(axis='both', which='minor', labelsize=16, width=1, size=2)
axs[-1].text(0.975, 0.1, "average telluric model", horizontalalignment='right', verticalalignment='bottom', transform=axs[-1].transAxes, fontsize=15)

axs[2].set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\mathring{A}^{-1}$)", fontsize=18)
axs[len(name_list)].yaxis.set_label_coords(-0.04, 0)

fig.subplots_adjust(hspace=0)
# plt.show()
fig.savefig('nires_1d2d.pdf')
