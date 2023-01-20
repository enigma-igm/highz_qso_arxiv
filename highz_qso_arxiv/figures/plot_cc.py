import os
import corner
import numpy as np
import matplotlib as mpl
import astropy.units as u
import mpl_scatter_density
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

from astropy.io import ascii
from astropy.table import Table
from highz_qso_arxiv import plot

from highz_qso_arxiv.plot import plot_contour, plot_cline
from highz_qso_arxiv.util import get_project_root, inverse
from highz_qso_arxiv.resource.filters import ukirt_J, wise_W1, decam_z

from IPython import embed

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--no-litqso', action='store_true')
parser.add_argument('--no-litdwf', action='store_true')
parser.add_argument('--no-colorcut', action='store_true')
parser.add_argument('--no-cov', action='store_true')
parser.add_argument('--no-meancov', action='store_true')
parser.add_argument('--show', action='store_true')
args = parser.parse_args()

def flux_to_mag(flux):
    return 22.5 - 2.5 * np.log10(flux)
def mag_to_flux(mag):
    return 10**((22.5 - mag) / 2.5)
def flux_ratio(color):
    return 10**(-color/2.5)

### Load Data

# 1. Riccardo's targets, contour plot

# all ukidss targets
uki_all = pyfits.getdata('../resource/catalog/UKIDSS_catalog_clean.fits', 1)
fz_uki, fz_err_uki = uki_all['flux_z'], uki_all['flux_err_z']
fW1_uki, fW1_err_uki = uki_all['flux_W1'], uki_all['flux_err_W1']
fJ_uki, fJ_err_uki = uki_all['flux_J'], uki_all['flux_J_err']

# TODO: add non-detection
mz_uki, mW1_uki, mJ_uki = flux_to_mag(fz_uki), flux_to_mag(fW1_uki), flux_to_mag(fJ_uki)

# all viking targets
vik_all = pyfits.getdata('../resource/catalog/VIKING_catalog_clean_nobright.fits', 1)
fz_vik, fz_err_vik = vik_all['flux_z'], vik_all['flux_z_err']
fW1_vik, fW1_err_vik = vik_all['flux_w1'], vik_all['flux_w1_err']
fJ_vik, fJ_err_vik = vik_all['J_flux_aper_3p0'], vik_all['J_flux_aper_err_3p0']

# TODO: add non-detection
mz_vik, mW1_vik, mJ_vik = flux_to_mag(fz_vik), flux_to_mag(fW1_vik), flux_to_mag(fJ_vik)

mz_all = np.concatenate([mz_uki, mz_vik])
mW1_all = np.concatenate([mW1_uki, mW1_vik])
mJ_all = np.concatenate([mJ_uki, mJ_vik])
mask = np.isfinite(mz_all) & np.isfinite(mW1_all) & np.isfinite(mJ_all)
fz_all = np.concatenate([fz_uki, fz_vik])
fW1_all = np.concatenate([fW1_uki, fW1_vik])
fJ_all = np.concatenate([fJ_uki, fJ_vik])

# 2. Literature dwarfs (provided by Feige)

dwarf_cat = pyfits.getdata('../resource/catalog/mltq_flux_lsdr9.fits', 1)
fz_mltq, fz_err_mltq = dwarf_cat['flux_z'], np.sqrt(inverse(dwarf_cat['flux_ivar_z']))
fW1_mltq, fW1_err_mltq = dwarf_cat['flux_w1'], np.sqrt(inverse(dwarf_cat['flux_ivar_w1']))
mJ_mltq = dwarf_cat['t1_japermag3'] + 0.938
mJ_err_mltq = dwarf_cat['t1_japermag3err']
# TODO: check error propagation formula for the error
fJ_mltq, fJ_err_mltq = mag_to_flux(mJ_mltq), mag_to_flux(mJ_err_mltq)

mz_mltq = flux_to_mag(dwarf_cat['flux_z'])
mW1_mltq = flux_to_mag(dwarf_cat['flux_w1'])
type = dwarf_cat['t1_redshift']
# indm = np.where((type == 'M') & (mz - mJ < 2.5)) # Feige did this but I don't know why
indm = np.where(type == 'M') 
indl = np.where(type == 'L')
indt = np.where(type == 'T')

# 3. Dwarfs in our sample

cand_uki_cat = pyfits.getdata('../resource/catalog/hizqa_uki.fits', 1)
fz_dwarf_uki, fz_err_dwarf_uki = cand_uki_cat['flux_z'], cand_uki_cat['flux_err_z'] 
fW1_dwarf_uki, fW1_err_dwarf_uki = cand_uki_cat['flux_W1'], cand_uki_cat['flux_err_W1']
fJ_dwarf_uki, fJ_err_dwarf_uki = cand_uki_cat['flux_J'], cand_uki_cat['flux_J_err']
mz_dwarf_uki, mW1_dwarf_uki, mJ_dwarf_uki = flux_to_mag(fz_dwarf_uki), flux_to_mag(fW1_dwarf_uki), flux_to_mag(fJ_dwarf_uki)
label_uki = cand_uki_cat['label']
type_uki = list(cand_uki_cat['spectral type'])
type_init_uki = np.array([x[:1] for x in type_uki])
source_uki = cand_uki_cat['source']

cand_vik_cat = pyfits.getdata('../resource/catalog/hizqa_vik.fits', 1)
fz_dwarf_vik, fz_err_dwarf_vik = cand_vik_cat['flux_z'], cand_vik_cat['flux_z_err']
fW1_dwarf_vik, fW1_err_dwarf_vik = cand_vik_cat['flux_w1'], cand_vik_cat['flux_w1_err']
fJ_dwarf_vik, fJ_err_dwarf_vik = cand_vik_cat['J_flux_aper_3p0'], cand_vik_cat['J_flux_aper_err_3p0']
mz_dwarf_vik, mW1_dwarf_vik, mJ_dwarf_vik = flux_to_mag(fz_dwarf_vik), flux_to_mag(fW1_dwarf_vik), flux_to_mag(fJ_dwarf_vik)
label_vik = cand_vik_cat['label']
type_vik = list(cand_vik_cat['spectral type'])
type_init_vik = np.array([x[:1] for x in type_vik])
source_vik = cand_vik_cat['source']

mz_dwarf = np.concatenate((mz_dwarf_uki, mz_dwarf_vik))
mW1_dwarf = np.concatenate((mW1_dwarf_uki, mW1_dwarf_vik))
mJ_dwarf = np.concatenate((mJ_dwarf_uki, mJ_dwarf_vik))
type_init = np.concatenate((type_init_uki, type_init_vik))
fz_dwarf = np.concatenate((fz_dwarf_uki, fz_dwarf_vik))
fW1_dwarf = np.concatenate((fW1_dwarf_uki, fW1_dwarf_vik))
fJ_dwarf = np.concatenate((fJ_dwarf_uki, fJ_dwarf_vik))
fz_err_dwarf = np.concatenate((fz_err_dwarf_uki, fz_err_dwarf_vik))
fW1_err_dwarf = np.concatenate((fW1_err_dwarf_uki, fW1_err_dwarf_vik))
fJ_err_dwarf = np.concatenate((fJ_err_dwarf_uki, fJ_err_dwarf_vik))
source_dwarf = np.concatenate((source_uki, source_vik))
label_dwarf = np.concatenate((label_uki, label_vik))

mask_Mdwarf = (type_init == 'M') 
mask_Ldwarf = (type_init == 'L')
mask_Tdwarf = (type_init == 'T')
mask_unknown = (label_dwarf == 'uNQ')
mask_qso = (label_dwarf == 'QSO') & (source_dwarf != 'known quasar')

# 4. Known quasars
uki_litqso = pyfits.getdata('../resource/catalog/qso_literature_uki.fits', 1)
redshift_litqso_uki = uki_litqso['redshift']
fz_litqso_uki, fz_err_litqso_uki = uki_litqso['flux_z_2'], uki_litqso['flux_err_z']
fW1_litqso_uki, fW1_err_litqso_uki = uki_litqso['flux_W1_2'], uki_litqso['flux_err_W1']
fJ_litqso_uki, fJ_err_litqso_uki = uki_litqso['flux_J'], uki_litqso['flux_J_err']
mz_litqso_uki, mW1_litqso_uki, mJ_litqso_uki = flux_to_mag(fz_litqso_uki), flux_to_mag(fW1_litqso_uki), flux_to_mag(fJ_litqso_uki)

vik_litqso = pyfits.getdata('../resource/catalog/qso_literature_vik.fits', 1)
redshift_litqso_vik = vik_litqso['redshift']
fz_litqso_vik, fz_err_litqso_vik = vik_litqso['flux_z_2'], vik_litqso['flux_z_err']
fW1_litqso_vik, fW1_err_litqso_vik = vik_litqso['flux_w1_2'], vik_litqso['flux_w1_err']
fJ_litqso_vik, fJ_err_litqso_vik = vik_litqso['J_flux_aper_3p0'], vik_litqso['J_flux_aper_err_3p0']
mz_litqso_vik, mW1_litqso_vik, mJ_litqso_vik = flux_to_mag(fz_litqso_vik), flux_to_mag(fW1_litqso_vik), flux_to_mag(fJ_litqso_uki)

redshift_litqso = np.concatenate((redshift_litqso_uki, redshift_litqso_vik))
mz_litqso = np.concatenate((mz_litqso_uki, mz_litqso_vik))
mW1_litqso = np.concatenate((mW1_litqso_uki, mW1_litqso_vik))
mJ_litqso = np.concatenate((mJ_litqso_uki, mJ_litqso_vik))
fz_litqso = np.concatenate((fz_litqso_uki, fz_litqso_vik))
fW1_litqso = np.concatenate((fW1_litqso_uki, fW1_litqso_vik))
fJ_litqso = np.concatenate((fJ_litqso_uki, fJ_litqso_vik))

# 5. Quasar redshift track
sim_qso = pyfits.getdata('../resource/catalog/high_z_QSOs.fits', 1)
redshift_simqso = sim_qso['z']
redshift_bins = np.linspace(6, 8, 20)
bins_center = (redshift_bins[1:] + redshift_bins[:-1]) / 2
bins_idx = np.digitize(redshift_simqso, redshift_bins)
mz_simqso = flux_to_mag(sim_qso['DECam-DECaLS-z'])
mW1_simqso = flux_to_mag(sim_qso['WISE-unWISE-W1'])
mJ_simqso = flux_to_mag(sim_qso['VISTA-VISTA-J'])
zJ_simqso = mz_simqso - mJ_simqso
JW1_simqso = mJ_simqso - mW1_simqso
# calculate the average zJ of each bin
zJ_simqso_mean = np.zeros(len(redshift_bins) - 1)
JW1_simqso_mean = np.zeros(len(redshift_bins) - 1)
for i in range(len(redshift_bins) - 1):
    zJ_simqso_mean[i] = np.mean(zJ_simqso[bins_idx == i + 1])
    JW1_simqso_mean[i] = np.mean(JW1_simqso[bins_idx == i + 1])

# 6. Dwarf spectral type track
czJ_dwarftrack, cJW1_dwarftrack = [], []
temp = []
setting = 'logg-5-metallicity-0'
btsettl_path = get_project_root() / 'resource' / 'dwarf' / 'bt-settl'
files = [f for f in os.listdir(btsettl_path / setting) if f.startswith('lte-')]
for f in files:
    dat = pyfits.getdata(btsettl_path / setting / f, 1)
    wl, flux = dat['wave'], dat['flux']
    wl, ind = np.unique(wl, return_index=True)
    flux = flux[ind]

    wl, flux = wl * u.AA, flux * u.erg / u.cm ** 2/ u.s / u.AA

    _zJ = decam_z.get_ab_magnitudes(flux, wl)[0][0] - ukirt_J.get_ab_magnitudes(flux, wl)[0][0]
    _JW1 = ukirt_J.get_ab_magnitudes(flux, wl)[0][0] - wise_W1.get_ab_magnitudes(flux, wl)[0][0]
    czJ_dwarftrack.append(_zJ)
    cJW1_dwarftrack.append(_JW1)
    # extract the number between 'lte-' and 'K.fits'
    temp.append(int(f[4:-8]))
index = np.argsort(temp)
czJ_dwarftrack = np.array(czJ_dwarftrack)[index]
cJW1_dwarftrack = np.array(cJW1_dwarftrack)[index]

dwarf_color = [["M5", 0.91, 0.47, 0.55, 0.45, 0.32, 0.11, 0.17],
                ["M6", 1.45, 0.60, 0.67, 0.53, 0.39, 0.22, 0.21],
                ["M7", 1.77, 0.70, 0.78, 0.56, 0.44, 0.25, 0.24],
                ["M8", 1.93, 0.77, 0.87, 0.58, 0.47, 0.26, 0.26],
                ["M9", 1.99, 0.82, 0.96, 0.60, 0.51, 0.27, 0.27],
                ["L0", 2.01, 0.86, 1.04, 0.63, 0.54, 0.29, 0.27],
                ["L1", 2.02, 0.88, 1.11, 0.67, 0.58, 0.33, 0.28],
                ["L2", 2.04, 0.90, 1.18, 0.73, 0.63, 0.40, 0.28],
                ["L3", 2.10, 0.92, 1.23, 0.79, 0.67, 0.48, 0.29],
                ["L4", 2.20, 0.94, 1.27, 0.86, 0.71, 0.56, 0.30],
                ["L5", 2.33, 0.97, 1.31, 0.91, 0.74, 0.65, 0.32],
                ["L6", 2.51, 1.00, 1.33, 0.96, 0.75, 0.72, 0.36],
                ["L7", 2.71, 1.04, 1.35, 0.97, 0.75, 0.77, 0.41],
                ["L8", 2.93, 1.09, 1.21, 0.96, 0.71, 0.79, 0.48],
                ["L9", 3.15, 1.16, 1.20, 0.90, 0.65, 0.79, 0.57],
                ["T0", 3.36, 1.23, 1.19, 0.80, 0.56, 0.76, 0.68],
                ["T1", 3.55, 1.33, 1.19, 0.65, 0.45, 0.71, 0.82],
                ["T2", 3.70, 1.43, 1.18, 0.46, 0.31, 0.65, 0.99],
                ["T3", 3.82, 1.55, 1.18, 0.25, 0.16, 0.59, 1.19],
                ["T4", 3.90, 1.68, 1.17, 0.02, 0.01, 0.55, 1.43],
                ["T5", 3.95, 1.81, 1.16, -0.19, -0.11, 0.54, 1.70],
                ["T6", 3.98, 1.96, 1.16, -0.35, -0.19, 0.59, 2.02],
                ["T7", 4.01, 2.11, 1.15, -0.43, -0.20, 0.70, 2.38],
                ["T8", 4.08, 2.26, 1.15, -0.36, -0.09, 0.90, 2.79]]
spectype = [dwarf_color[i][0] for i in range(len(dwarf_color))]
czY = [dwarf_color[i][2] for i in range(len(dwarf_color))]
cYJ = [dwarf_color[i][3] for i in range(len(dwarf_color))]
cJH = [dwarf_color[i][4] for i in range(len(dwarf_color))]
cHK = [dwarf_color[i][5] for i in range(len(dwarf_color))]
cKW1 = [dwarf_color[i][6] for i in range(len(dwarf_color))]
cW1W2 = [dwarf_color[i][7] for i in range(len(dwarf_color))]
czJ = np.array(czY) + np.array(cYJ) + 0.54 - 0.938
cW1J = -np.array(cKW1) - np.array(cHK) - np.array(cJH) + 2.699 - 0.938
cJW1 = -cW1J

fig, ax = plt.subplots(1, 1, figsize=(10, 10))

ii = 0

# plotting Riccardo's targets
ax = plot_contour(mz_all[mask] - mJ_all[mask],
                  mJ_all[mask] - mW1_all[mask], 
                  ax=ax, bins=150, 
                  levels=np.linspace(0.1, 0.9, 5),
                  color='grey', zorder=0)
ax = plot_contour(mz_simqso - mJ_simqso,
                  mJ_simqso - mW1_simqso,
                  ax=ax, bins=50,
                  levels=np.linspace(0.1, 0.9, 5),
                  color='tomato', zorder=1)
# plotting literature dwarfs
ax.scatter(mz_mltq[indm] - mJ_mltq[indm], mJ_mltq[indm] - mW1_mltq[indm], zorder=2,
           s=20, marker='o', facecolors='none', edgecolors=CB_color_cycle[ii], alpha=0.6, label=r'$\textbf{lit. M dwarfs}$')
ax.scatter(mz_mltq[indl] - mJ_mltq[indl], mJ_mltq[indl] - mW1_mltq[indl], zorder=2,
           s=20, marker='v', facecolors='none', edgecolors=CB_color_cycle[ii], alpha=0.6, label=r'$\textbf{lit. L dwarfs}$')
ax.scatter(mz_mltq[indt] - mJ_mltq[indt], mJ_mltq[indt] - mW1_mltq[indt], zorder=2,
           s=20, marker='s', facecolors='none', edgecolors=CB_color_cycle[ii], alpha=0.6, label=r'$\textbf{lit. T dwarfs}$')
ax.text(-0.5, -1.4, r'$\textbf{M dwarf}$', fontsize=15, zorder=2, color=CB_color_cycle[ii])
ax.text(2.2, 1., r'$\textbf{M dwarf}$', fontsize=15, zorder=2, color=CB_color_cycle[ii])
ax.text(4, -2, r'$\textbf{M dwarf}$', fontsize=15, zorder=2, color=CB_color_cycle[ii])

# plotting our dwarfs
ax.scatter(mz_dwarf[mask_Mdwarf] - mJ_dwarf[mask_Mdwarf], 
           mJ_dwarf[mask_Mdwarf] - mW1_dwarf[mask_Mdwarf], 
           s=30, facecolors=CB_color_cycle[ii], marker='o', 
           edgecolors='black',
           label=r'$\textbf{M dwarf}$', zorder=6)
ax.scatter(mz_dwarf[mask_Ldwarf] - mJ_dwarf[mask_Ldwarf], 
           mJ_dwarf[mask_Ldwarf] - mW1_dwarf[mask_Ldwarf], 
           s=30, facecolors=CB_color_cycle[ii], marker='v', 
           edgecolors='black',
           label=r'$\textbf{L dwarf}$', zorder=6)
ax.scatter(mz_dwarf[mask_Tdwarf] - mJ_dwarf[mask_Tdwarf],
           mJ_dwarf[mask_Tdwarf] - mW1_dwarf[mask_Tdwarf],
           s=30, facecolors=CB_color_cycle[ii], marker='s', 
           edgecolors='black',
           label=r'$\textbf{T dwarf}$', zorder=6)
ax.scatter(mz_dwarf[mask_unknown] - mJ_dwarf[mask_unknown], 
           mJ_dwarf[mask_unknown] - mW1_dwarf[mask_unknown], 
           s=20, c='black', marker='d', label=r'$\textbf{unknown non-QSO}$', 
           zorder=5, alpha=0.8)
ax.scatter(mz_dwarf[mask_qso] - mJ_dwarf[mask_qso], 
           mJ_dwarf[mask_qso] - mW1_dwarf[mask_qso], 
           s=100, c='red', marker='*', label=r'$\textbf{QSO}$', zorder=5)

# plotting literature quasars
ax.scatter(mz_litqso - mJ_litqso, mJ_litqso - mW1_litqso, s=40, 
           facecolors='none', edgecolors='red', alpha=0.5, 
           marker='*', label=r'$\textbf{Known quasars}$', zorder=4)
z7mask = redshift_litqso > 7.
ax.scatter(mz_litqso[z7mask] - mJ_litqso[z7mask], mJ_litqso[z7mask] - mW1_litqso[z7mask], s=100,
           facecolors='none', edgecolors='red', alpha=1.0, 
           marker='*', label=r'$\textbf{Known quasars}$ ${\bf (z>7)}$', zorder=4)

# plotting quasar track
plot_cline(zJ_simqso_mean, JW1_simqso_mean, ax=ax, cmap=plt.get_cmap('autumn_r'), linewidth=2)

# plotting dwarf track
# ax.scatter(czJ_dwarftrack, cJW1_dwarftrack, c=CB_color_cycle[2], s=5, alpha=0.5, zorder=2)
ax.plot(czJ_dwarftrack, cJW1_dwarftrack, alpha=0.8, 
        c=CB_color_cycle[ii], zorder=2, label=r'$\textbf{Dwarf track (Allard+2012)}$')
ax.plot(czJ, cJW1, c='black', alpha=0.6, zorder=2, label=r'$\textbf{Dwarf track (Skrzypek+2015)}$')

# plotting Feige's color-color cuts
ax.text(-0.8, 1.5, r'${\bf z\lesssim 6.5}$', fontsize=12, color='black', zorder=5)
ax.text(0.8, 1.5, r'${\bf z\gtrsim 6.5}$', fontsize=12, color='black', zorder=5)
ax.text(7, -0.6, r'${\bf z>7}$', fontsize=12, color='black', zorder=5)

ax.plot([-1., 0.562], [-0.261, -0.261], ls='dashed', lw=2, alpha=0.6, color='blue', zorder=2)
ax.plot([-1., 0.562], [1.239, 1.239], ls='dashed', lw=2, alpha=0.6, color='blue', zorder=2)
ax.plot([0.562, 0.562], [-0.261, 1.239], ls='dashed', lw=2, alpha=0.6, color='blue', zorder=2)
ax.plot([0.562, 1.562], [1.239, 1.239], ls='dashed', lw=2, alpha=0.6, color='orange', zorder=2)
ax.plot([0.562, 1.062], [-0.261, -0.261], ls='dashed', lw=2, alpha=0.6, color='orange', zorder=2)
ax.plot([1.062, 1.562], [-0.261, 1.239], ls='dashed', lw=2, alpha=0.6, color='orange', zorder=2)
ax.plot([2.062, 2.062], [-0.261, 1.239], ls='dashed', lw=2, alpha=0.6, color='red', zorder=2)
ax.plot([2.062, 8], [-0.261, -0.261], ls='dashed', lw=2, alpha=0.6, color='red', zorder=2)
ax.plot([2.062, 8], [1.239, 1.239], ls='dashed', lw=2, alpha=0.6, color='red', zorder=2)

ax.legend(loc='upper right', fontsize=12, frameon=True, framealpha=0.8)
ax.set_xlabel(r'$z_{\rm AB}$ - $J_{\rm AB}$', fontsize=25)
ax.set_ylabel(r'$J_{\rm AB}$ - $W1_{\rm AB}$', fontsize=25)
ax.set_xlim(-1, 8)
ax.set_ylim(-4, 6)

ax.tick_params(axis='both', which='major', labelsize=20)

if args.show:
    plt.show()
else:
    plt.savefig('color_color.png', dpi=300)