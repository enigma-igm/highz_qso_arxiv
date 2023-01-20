import os
import speclite
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

from astropy.io import fits, ascii
from highz_qso_arxiv.util import get_project_root
from highz_qso_arxiv.resource.filters import ukirt_J

from IPython import embed

import matplotlib as mpl
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
CB_color_cycle = ['#2166ac', '#b2182b']
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 

path = get_project_root()
ARXIV_PATH = os.path.join(path, "arxiv")
RESOURCE_PATH = os.path.join(path, "resource")

m_J = 21.

redshift = 7.
qso_spec = ascii.read(os.path.join(RESOURCE_PATH, "Selsing2015.dat"))
wl_qso_rest = qso_spec["col1"]
wl_qso_obs = wl_qso_rest * (1 + redshift)
flux_qso = qso_spec["col2"] # in 1e-17 erg/s/cm2/A
flux_qso_err = qso_spec["col3"]
wl_lya = 1215.67 * (1 + redshift)
trough = wl_qso_obs < wl_lya
flux_qso[trough] = 0.0
label_qso = rf"$z={redshift}$; Selsing+2016"

m_J_temp = ukirt_J.get_ab_magnitudes(flux_qso * 1e-17 * u.erg / u.s / u.cm**2 / u.AA, wl_qso_obs*u.AA)[0][0]
scale = 10**(-(m_J-m_J_temp)/2.5)
flux_qso = flux_qso * scale

star_type = 'M8'
star_spec = ascii.read(os.path.join(RESOURCE_PATH, f"dwarf/irtf_{star_type}_1.dat"))
wl_star_obs = star_spec["col1"] * 1e4
flux_star = star_spec["col2"] * u.W / u.m**2 / u.um # in 1e-17 erg/s/cm2/A
flux_star = flux_star.to(1e-17 * u.erg / u.s / u.cm**2 / u.AA).value
label_star = star_type
wl_star_obs = wl_star_obs[flux_star > 0]
flux_star = flux_star[flux_star > 0]

# flux, wl = ukirt_J.pad_spectrum(flux_star, wl_star_obs, method='zero')
m_J_temp = ukirt_J.get_ab_magnitudes(flux_star * 1e-17 * u.erg / u.s / u.cm**2 / u.AA, wl_star_obs*u.AA)[0][0]
scale = 10**(-(m_J-m_J_temp)/2.5)
flux_star = flux_star * scale

fig, ax = plt.subplots(figsize=(12,6))
# ax.plot(wl_qso_obs[~trough], flux_qso[~trough], label=label_qso, color="black", zorder=3)
# ax.plot(wl_qso_obs, flux_qso, color="grey", ls="--", zorder=2)
ax.plot(wl_qso_obs, flux_qso, label=label_qso, color=CB_color_cycle[1], alpha=0.8, lw=2, zorder=3)
# ax.plot(wl_qso_obs, flux_qso, color=CB_color_cycle[1], alpha=0.5, lw=2, ls="--", zorder=2)

ax.plot(wl_star_obs, flux_star, label=label_star, color=CB_color_cycle[0], alpha=0.8, zorder=1, lw=2)

# ax.fill_between(wl_obs, flux-flux_err, flux+flux_err, color="black", alpha=0.2)
# ax.set_xlim(7300, 10500)
ax.legend(loc="upper left", fontsize=20)
# ax.set_xlabel(r"wavelength ($\mathring{A}$)", fontsize=28)
ax.set_xlabel(r"wavelength ($\mu m$)", fontsize=28)
ax.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=28)

filter_J = ascii.read(os.path.join(RESOURCE_PATH, f"filter/UKIRT_UKIDSS.J.dat"))
filter_z = ascii.read(os.path.join(RESOURCE_PATH, f"filter/decam_z.dat"))
wave_J = filter_J["col1"]
tran_J = filter_J["col2"]
wave_z = filter_z["col1"]
tran_z = filter_z["col2"]
# ax.plot(wave_J, tran_J, color="blue", alpha=0.2, zorder=0)
# ax.fill_between(wave_J, tran_J, color="blue", alpha=0.2, zorder=0)
# ax.text(np.mean(wave_J), 0.05, r"$J$", fontsize=25, color="black")
# ax.plot(wave_z, tran_z, color="red", alpha=0.2, zorder=0)
# ax.fill_between(wave_z, tran_z, color="red", alpha=0.2, zorder=0)
# ax.text(np.mean(wave_z), 0.06, r"$z$", fontsize=25, color="black")

ax.plot(wave_J, tran_J*2, color="black", alpha=0.2, zorder=0)
ax.text(np.mean(wave_J), 1.2, r"$J$", fontsize=25, color="black")
ax.plot(wave_z, tran_z*2, color="black", alpha=0.2, zorder=0)
ax.text(np.mean(wave_z), 1.2, r"$z$", fontsize=25, color="black")

ax.set_xlim(8000, 14000)
ax.set_ylim(0, 1.5)

# put text on upper right corner
ax.text(0.95, 0.95, r"$m_{J}=21.0$", ha="right", va="top", transform=ax.transAxes, fontsize=25)
ax.tick_params(labelsize=20)
ax.set_ylim(-0.0, 1.8)
ax.set_yticks([0.5, 1.0, 1.5])
ax.tick_params(axis='both', which='major', labelsize=30, width=1, size=6)
ax.tick_params(axis='both', which='minor', labelsize=30, width=1, size=3)
ax.set_xticklabels(["0.8", "0.9", "1.0", "1.1", "1.2", "1.3", "1.4"])

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--show', action='store_true')
args = parser.parse_args()

if args.show:
    plt.show()
else:
    plt.savefig('qso_star.pdf')