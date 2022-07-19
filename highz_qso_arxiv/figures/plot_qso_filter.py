import numpy as np
import matplotlib as mpl
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from matplotlib.colors import LinearSegmentedColormap

from highz_qso_arxiv import plot as hz_plot

from IPython import embed

fig, ax = plt.subplots(figsize=(16, 8))

color = plt.cm.bwr(np.linspace(0.2,0.8,7))
colors = [color[0], color[-1]] # first color is black, last is red
cm = LinearSegmentedColormap.from_list("Custom", colors, N=7)
color = cm(np.linspace(0,1,7))

dat = ascii.read("../resource/filter/decam_z.dat")
wave_decam_z = dat["Lambda"]
transmission_decam_z = dat["Transmission"]

dat = ascii.read("../resource/filter/decam_r.dat")
wave_decam_r = dat["Lambda"]
transmission_decam_r = dat["Transmission"]

dat = ascii.read("../resource/filter/vista_J.dat")
wave_vista_J = dat["Lambda"]
transmission_vista_J = dat["Transmission"]

for i, redshift in enumerate(np.arange(5.5, 8.5, 0.5)):
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
    ax.fill_between(wl_obs, 0, flux, color=color[i], zorder=10-i, label=label, alpha=1)

    # ymin, ymax = ax.get_ylim()
ax.set_ylim(10, 120)
ax.set_yticklabels([])
ax.set_xlim(wave_decam_r[0], wave_vista_J[-1])
ax.legend(loc="center left", fontsize=15)
ax.set_xlabel("Wavelength [Angstrom]", fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)

color = cm(np.linspace(0,1,3))
ax2 = ax.twinx()
ax2.set_ylim(3, 0)
ax2.set_yticklabels([])
ax2.plot(wave_decam_r, transmission_decam_r, color=color[0], label="DECam r")
ax2.plot(wave_decam_z, transmission_decam_z, color=color[1], label="DECam z")
ax2.plot(wave_vista_J, transmission_vista_J, color=color[2], label="VISTA J")
ax2.legend(fontsize=15)

# plt.show()
plt.savefig("./qso_filter.pdf")