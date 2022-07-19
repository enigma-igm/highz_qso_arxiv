import numpy as np
import matplotlib as mpl
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits, ascii
from matplotlib.colors import LinearSegmentedColormap

from highz_qso_arxiv import plot as hz_plot
from highz_qso_arxiv.util import inverse, ivarsmooth

from IPython import embed

# fig, ax = plt.subplots(figsize=(20, 8))

# color = plt.cm.coolwarm(np.linspace(0.05,0.95,7))
# colors = [color[0], color[-1]] # first color is black, last is red
# cm = LinearSegmentedColormap.from_list("Custom", colors, N=7)
# color = cm(np.linspace(0,1,7))

# dat = ascii.read("../resource/filter/decam_z.dat")
# wave_decam_z = dat["Lambda"]
# transmission_decam_z = dat["Transmission"]

# dat = ascii.read("../resource/filter/decam_r.dat")
# wave_decam_r = dat["Lambda"]
# transmission_decam_r = dat["Transmission"]

# # dat = ascii.read("../resource/filter/vista_J.dat")
# # wave_vista_J = dat["Lambda"]
# # transmission_vista_J = dat["Transmission"]

# fits_file="../arxiv/LRIS_2203/LRIS_220306/reduced/all/J0901+2906/J0901+2906_coadd_tellcorr.fits"
# hdul = fits.open(fits_file)
# data = Table(hdul[1].data)
# try:
#     wave, flux, flux_ivar = data["wave"], data["flux"], data["ivar"]
# except KeyError:
#     wave, flux, flux_ivar = data["OPT_WAVE"], data["OPT_FLAM"], data["OPT_FLAM_IVAR"]

# mask = wave>1.
# flux_sm, flux_ivar_sm = ivarsmooth(flux, flux_ivar, 9)
# flux_sm = flux_sm / np.mean(flux_sm)
# ax.plot(wave[mask], flux_sm[mask], label="z=6.1", color="black", lw=1.5, zorder=5)

# for dwarf in ["M8", "L5"]:
#     dat = ascii.read(f"../resource/dwarf/keck_lris_{dwarf}_1.dat")
#     wave = dat["col1"]
#     flux = dat["col2"]
#     flux = flux / np.mean(flux) / 2
#     ax.plot(wave, flux, label=dwarf)
# ax.set_xlim(5400, 10200)
# ax.set_ylim(-0.1, 4)
# ax.set_yticklabels([])
# ax.legend()

# color = cm(np.linspace(0,1,2))
# ax2 = ax.twinx()
# ax2.set_ylim(3, 0)
# ax2.set_yticklabels([])
# ax2.plot(wave_decam_r, transmission_decam_r, color=color[0], label="DECam r", zorder=1)
# ax2.plot(wave_decam_z, transmission_decam_z, color=color[1], label="DECam z", zorder=1)
# # ax2.plot(wave_vista_J, transmission_vista_J, color=color[2], label="VISTA J")
# ax2.legend(fontsize=15)

# plt.show()


fig, ax = plt.subplots(figsize=(16, 8))

color = plt.cm.bwr(np.linspace(0.2,0.8,4))
colors = [color[0], color[-1]] # first color is black, last is red
cm = LinearSegmentedColormap.from_list("Custom", colors, N=4)
color = cm(np.linspace(0,1,4))

dat = ascii.read("../resource/filter/decam_z.dat")
wave_decam_z = dat["Lambda"]
transmission_decam_z = dat["Transmission"]

dat = ascii.read("../resource/filter/decam_r.dat")
wave_decam_r = dat["Lambda"]
transmission_decam_r = dat["Transmission"]

dat = ascii.read("../resource/filter/vista_J.dat")
wave_vista_J = dat["Lambda"]
transmission_vista_J = dat["Transmission"]

for i, redshift in enumerate(np.arange(5.5, 7.5, 0.5)):
    dat = ascii.read("../resource/Selsing2015.dat")
    wl_rest = dat["col1"]
    wl_obs = wl_rest * (1 + redshift)
    flux = dat["col2"] # in 1e-17 erg/s/cm2/A
    wl_lya = 1215.67 * (1 + redshift)
    trough = wl_obs < wl_lya

    flux = flux + i*10
    flux[trough] = 0
    label = rf"$z={redshift}$"

    ax.plot(wl_obs, flux, color=color[i], zorder=10-i)
    ax.fill_between(wl_obs, 0, flux, color=color[i], zorder=10-i, label=label)

    # ymin, ymax = ax.get_ylim()
ax.set_ylim(10, 120)
ax.set_yticklabels([])
ax.set_xlim(wave_decam_r[0], wave_decam_z[-1])
ax.set_xlabel("Wavelength [Angstrom]", fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)

# styles = ["-", "dotted"]
alphas = [0.5, 1.0]
colors = ["black", "grey"]
for i, dwarf in enumerate(["M8", "L5"]):
    dat = ascii.read(f"../resource/dwarf/keck_lris_{dwarf}_1.dat")
    wave = dat["col1"]
    flux = dat["col2"]
    flux = 15 * flux / np.mean(flux) + 10
    ax.plot(wave, flux, label=dwarf, lw=2, zorder=20, color="black", alpha=alphas[i])

ax.legend(loc="center left", fontsize=15)

color = cm(np.linspace(0,1,2))
ax2 = ax.twinx()
ax2.set_ylim(3, 0)
ax2.set_yticklabels([])
ax2.plot(wave_decam_r, transmission_decam_r, color=color[0], label="DECam r")
ax2.plot(wave_decam_z, transmission_decam_z, color=color[1], label="DECam z")
# ax2.plot(wave_vista_J, transmission_vista_J, color=color[2], label="VISTA J")
ax2.legend(fontsize=15)

# plt.show()
plt.savefig("contaminant.pdf")