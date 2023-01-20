import os
import corner
import argparse
import numpy as np
import matplotlib as mpl
import astropy.units as u
import mpl_scatter_density
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

from astropy.io import ascii
from astropy.table import Table
from highz_qso_arxiv import plot

from highz_qso_arxiv.util.photutil import add_noise
from highz_qso_arxiv.util import get_project_root, inverse
from highz_qso_arxiv.plot import plot_contour, plot_cov, plot_cline, error_cov
from highz_qso_arxiv.resource.filters import ukirt_J, wise_W1, decam_z

from IPython import embed

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 

def flux_to_mag(flux):
    return 22.5 - 2.5 * np.log10(flux)
def mag_to_flux(mag):
    return 10**((22.5 - mag) / 2.5)
def flux_ratio(color):
    return 10**(-color/2.5)
def flux_ratio_cov(X, Y, Z, Xerr, Yerr, Zerr):
    # X/Z vs. Y/Z
    cov = np.zeros((len(X), 2, 2))
    cov[:, 0, 0] = (X / Z)**2 * ((Xerr / X)**2 + (Zerr / Z)**2)
    cov[:, 1, 1] = (Y / Z)**2 * ((Yerr / Y)**2 + (Zerr / Z)**2)
    cov[:, 0, 1] = cov[:, 1, 0] = X * Y * Zerr**2 / Z**4
    return cov

# Dwarf spectral type track
def dwarf_track(setting):
    czJ_dwarftrack, cJW1_dwarftrack = [], []
    temp = []
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
    return temp, czJ_dwarftrack, cJW1_dwarftrack
temp, czJ_dwarftrack, cJW1_dwarftrack = dwarf_track('logg-5-metallicity-0')

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
cW2J = -np.array(cW1W2) - np.array(cKW1) - np.array(cHK) - np.array(cJH) + 3.339 - 0.938
cYJ = np.array(cYJ) + 0.634 - 0.938
cHJ = -np.array(cJH) + 1.379 - 0.938
cKJ = -np.array(cHK) - np.array(cJH) + 1.900 - 0.938
cW1W2 = np.array(cW1W2) + 2.699 - 3.339
fzJ = flux_ratio(czJ)
fW1J = flux_ratio(cW1J)
fW2J = flux_ratio(cW1W2)
fYJ = flux_ratio(cYJ)
fHJ = flux_ratio(cHJ)
fKJ = flux_ratio(cKJ)

from tqdm import trange
all_syn_flux = np.array([], dtype=np.int64).reshape(7,0)
dummy_fJ = mag_to_flux(21)
syn_flux = np.zeros((7, len(fzJ)))
syn_flux[0] = dummy_fJ
syn_flux[1] = fzJ * dummy_fJ
syn_flux[2] = fYJ * dummy_fJ
syn_flux[3] = fHJ * dummy_fJ
syn_flux[4] = fKJ * dummy_fJ
syn_flux[5] = fW1J * dummy_fJ
syn_flux[6] = fW2J * dummy_fJ
grid_zJ = np.linspace(np.min(syn_flux[1]), np.max(syn_flux[1]), 100000)

# assign each point a spectral type
syn_spectype = []
for i in range(len(grid_zJ)):
    idx = np.argmin(np.abs(grid_zJ[i] - syn_flux[1]))
    syn_spectype.append(spectype[idx])
syn_spectype = np.array(syn_spectype)

_syn_flux = np.zeros((7, len(grid_zJ)))
_syn_flux[0] = dummy_fJ
_syn_flux[1] = grid_zJ
# interpolate the others
from scipy.interpolate import interp1d
for i in range(2, 7):
    f = interp1d(syn_flux[1], syn_flux[i], kind='cubic')
    _syn_flux[i] = f(grid_zJ)

syn_flux, syn_flux_err = add_noise(_syn_flux)
all_syn_flux = np.hstack((all_syn_flux, syn_flux))
np.save('../resource/catalog/all_syn_flux_noisy.npy', all_syn_flux)
np.save('../resource/catalog/all_syn_flux_err.npy', syn_flux_err)
np.save('../resource/catalog/all_syn_spectype.npy', syn_spectype)