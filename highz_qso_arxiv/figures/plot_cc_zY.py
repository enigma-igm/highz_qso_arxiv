import os
from re import A
import corner
import argparse
from matplotlib import markers
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

# CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
#                   '#f781bf', '#a65628', '#984ea3',
#                   '#999999', '#e41a1c', '#dede00']
CB_color_cycle = ['#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e', '#f781bf']

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 

parser = argparse.ArgumentParser()
parser.add_argument('--show', action='store_true')
parser.add_argument('--save', action='store_true')
args = parser.parse_args()

def flux_to_mag(flux, fluxerr, sigma=1):
    mask = flux > fluxerr * sigma
    mag = np.zeros_like(flux)
    mag[mask] = -2.5 * np.log10(flux[mask]) + 22.5
    mag[~mask] = -2.5 * np.log10(fluxerr[~mask] * sigma) + 22.5
    return mag
def fluxerr_to_magerr(fluxerr, flux, sigma=1):
    mask = flux > fluxerr * sigma
    magerr = np.zeros_like(fluxerr)
    magerr[mask] = 2.5 / np.log(10) * fluxerr[mask] / flux[mask]
    magerr[~mask] = 2.5 / np.log(10) * fluxerr[~mask] / (fluxerr[~mask] * sigma)
    return magerr
def mag_to_flux(mag):
    flux = np.zeros_like(mag)
    flux = 10**(-0.4 * (mag - 22.5))
    return flux
def magerr_to_fluxerr(magerr, flux):
    fluxerr = np.zeros_like(magerr)
    fluxerr = flux * np.log(10) / 2.5 * magerr
    return fluxerr
def flux_ratio(color):
    return 10**(-color/2.5)
def color_cov(X, Y, Z, Xerr, Yerr, Zerr):
    # X/Z vs. Y/Z
    cov = np.zeros((len(X), 2, 2))
    cov[:, 0, 0] = Xerr**2 + Zerr**2
    cov[:, 1, 1] = Yerr**2 + Zerr**2
    cov[:, 0, 1] = cov[:, 1, 0] = Zerr**2
    return cov

### Load Data

# 1. Riccardo's targets, contour plot

# all ukidss targets
uki_all = pyfits.getdata('../resource/catalog/UKIDSS_catalog_clean.fits', 1)
fz_uki, fz_err_uki = uki_all['flux_z'], uki_all['flux_err_z']
fW1_uki, fW1_err_uki = uki_all['flux_W1'], uki_all['flux_err_W1']
fW2_uki, fW2_err_uki = uki_all['flux_W2'], uki_all['flux_err_W2']
fY_uki, fY_err_uki = uki_all['Y_flux_aper_3p0'], uki_all['Y_flux_aper_err_3p0']
fJ_uki, fJ_err_uki = uki_all['flux_J'], uki_all['flux_J_err']
fH_uki, fH_err_uki = uki_all['H_flux_aper_3p0'], uki_all['H_flux_aper_err_3p0']
fK_uki, fK_err_uki = uki_all['K_flux_aper_3p0'], uki_all['K_flux_aper_err_3p0']

mz_uki = flux_to_mag(fz_uki, fz_err_uki)
mW1_uki = flux_to_mag(fW1_uki, fW1_err_uki)
mW2_uki = flux_to_mag(fW2_uki, fW2_err_uki)
mY_uki = flux_to_mag(fY_uki, fY_err_uki)
mJ_uki = flux_to_mag(fJ_uki, fJ_err_uki)
mH_uki = flux_to_mag(fH_uki, fH_err_uki)
mK_uki = flux_to_mag(fK_uki, fK_err_uki)

# all viking targets
vik_all = pyfits.getdata('../resource/catalog/VIKING_catalog_clean_nobright.fits', 1)
fz_vik, fz_err_vik = vik_all['flux_z'], vik_all['flux_z_err']
fW1_vik, fW1_err_vik = vik_all['flux_w1'], vik_all['flux_w1_err']
fW2_vik, fW2_err_vik = vik_all['flux_w2'], vik_all['flux_w2_err']
fY_vik, fY_err_vik = vik_all['Y_flux_aper_3p0'], vik_all['Y_flux_aper_err_3p0']
fJ_vik, fJ_err_vik = vik_all['J_flux_aper_3p0'], vik_all['J_flux_aper_err_3p0']
fH_vik, fH_err_vik = vik_all['H_flux_aper_3p0'], vik_all['H_flux_aper_err_3p0']
fK_vik, fK_err_vik = vik_all['K_flux_aper_3p0'], vik_all['K_flux_aper_err_3p0']

mz_vik = flux_to_mag(fz_vik, fz_err_vik)
mW1_vik = flux_to_mag(fW1_vik, fW1_err_vik)
mW2_vik = flux_to_mag(fW2_vik, fW2_err_vik)
mY_vik = flux_to_mag(fY_vik, fY_err_vik)
mJ_vik = flux_to_mag(fJ_vik, fJ_err_vik)
mH_vik = flux_to_mag(fH_vik, fH_err_vik)
mK_vik = flux_to_mag(fK_vik, fK_err_vik)

mz_all = np.concatenate([mz_uki, mz_vik])
mW1_all = np.concatenate([mW1_uki, mW1_vik])
mW2_all = np.concatenate([mW2_uki, mW2_vik])
mY_all = np.concatenate([mY_uki, mY_vik])
mJ_all = np.concatenate([mJ_uki, mJ_vik])
mH_all = np.concatenate([mH_uki, mH_vik])
mK_all = np.concatenate([mK_uki, mK_vik])
fz_all = np.concatenate([fz_uki, fz_vik])
fW1_all = np.concatenate([fW1_uki, fW1_vik])
fW2_all = np.concatenate([fW2_uki, fW2_vik])
fJ_all = np.concatenate([fJ_uki, fJ_vik])
fY_all = np.concatenate([fY_uki, fY_vik])
fH_all = np.concatenate([fH_uki, fH_vik])
fK_all = np.concatenate([fK_uki, fK_vik])

# 2. Literature dwarfs (provided by Feige)

# dwarf_cat = pyfits.getdata('../resource/catalog/mltq_flux_lsdr9.fits', 1)
# fz_mltq, fz_err_mltq = dwarf_cat['flux_z'], np.sqrt(inverse(dwarf_cat['flux_ivar_z']))
# fW1_mltq, fW1_err_mltq = dwarf_cat['flux_w1'], np.sqrt(inverse(dwarf_cat['flux_ivar_w1']))
# mJ_mltq = dwarf_cat['t1_japermag3'] + 0.938
# mJ_err_mltq = dwarf_cat['t1_japermag3err']
# # TODO: check error propagation formula for the error
# fJ_mltq, fJ_err_mltq = mag_to_flux(mJ_mltq), mag_to_flux(mJ_err_mltq)

# mz_mltq = flux_to_mag(dwarf_cat['flux_z'])
# mW1_mltq = flux_to_mag(dwarf_cat['flux_w1'])
# type = dwarf_cat['t1_redshift']
# # indm = np.where((type == 'M') & (mz - mJ < 2.5)) # Feige did this but I don't know why
# indm = np.where(type == 'M') 
# indl = np.where(type == 'L')
# indt = np.where(type == 'T')

# 3. Dwarfs in our sample

from astropy.io import fits
# cand_cat = Table(fits.getdata('../arxiv/hizqa_wildhunt.fits', 1))
cand_cat = Table(fits.getdata('../arxiv/hizqa_ric.fits', 1))

fz_cand, fz_err_cand = cand_cat['flux_z'], cand_cat['flux_err_z']
fW1_cand, fW1_err_cand = cand_cat['flux_W1'], cand_cat['flux_err_W1']
fW2_cand, fW2_err_cand = cand_cat['flux_W2'], cand_cat['flux_err_W2']
fJ_cand, fJ_err_cand = cand_cat['flux_J'], cand_cat['flux_err_J']
fY_cand, fY_err_cand = cand_cat['flux_Y'], cand_cat['flux_err_Y']
fH_cand, fH_err_cand = cand_cat['flux_H'], cand_cat['flux_err_H']
fK_cand, fK_err_cand = cand_cat['flux_K'], cand_cat['flux_err_K']

mz_cand = flux_to_mag(fz_cand, fz_err_cand)
mW1_cand = flux_to_mag(fW1_cand, fW1_err_cand)
mW2_cand = flux_to_mag(fW2_cand, fW2_err_cand)
mJ_cand = flux_to_mag(fJ_cand, fJ_err_cand)
mY_cand = flux_to_mag(fY_cand, fY_err_cand)
mH_cand = flux_to_mag(fH_cand, fH_err_cand)
mK_cand = flux_to_mag(fK_cand, fK_err_cand)

mz_err_cand = fluxerr_to_magerr(fz_err_cand, fz_cand)
mW1_err_cand = fluxerr_to_magerr(fW1_err_cand, fW1_cand)
mW2_err_cand = fluxerr_to_magerr(fW2_err_cand, fW2_cand)
mJ_err_cand = fluxerr_to_magerr(fJ_err_cand, fJ_cand)
mY_err_cand = fluxerr_to_magerr(fY_err_cand, fY_cand)
mH_err_cand = fluxerr_to_magerr(fH_err_cand, fH_cand)
mK_err_cand = fluxerr_to_magerr(fK_err_cand, fK_cand)

label_cand = cand_cat['label'].value
spectype_cand = cand_cat['spectral type'].value
run_cand = cand_cat['run'].value
type_init = np.array([x[:1] for x in spectype_cand])
instrument = np.array([x.split('-')[0] for x in run_cand])
source_cand = cand_cat['source']

mask_detected = (fz_cand > 0) & (fJ_cand > 0)
mask_Mdwarf = (type_init == 'M') & mask_detected
mask_Ldwarf = (type_init == 'L') & mask_detected
mask_Tdwarf = (type_init == 'T') & mask_detected
mask_dwarf = mask_Mdwarf | mask_Ldwarf | mask_Tdwarf
mask_unknown = (label_cand == 'uNQ') & mask_detected
mask_qso = (label_cand == 'QSO') & (source_cand != 'known quasar') & mask_detected

mask_mosfire = (instrument == 'mosfire') & mask_detected
mask_ric = (source_cand == 'RN') & mask_detected

cov_zY = color_cov(mz_cand[mask_detected], mY_cand[mask_detected], 
                    mJ_cand[mask_detected], mz_err_cand[mask_detected],
                    mY_err_cand[mask_detected], mJ_err_cand[mask_detected])
cov_dwf_zY = color_cov(mz_cand[mask_dwarf], mY_cand[mask_dwarf],
                        mJ_cand[mask_dwarf], mz_err_cand[mask_dwarf], 
                        mY_err_cand[mask_dwarf], mJ_err_cand[mask_dwarf])
cov_unk_zY = color_cov(mz_cand[mask_unknown], mY_cand[mask_unknown],
                        mJ_cand[mask_unknown], mz_err_cand[mask_unknown],
                        mY_err_cand[mask_unknown], mJ_err_cand[mask_unknown])
cov_qso_zY = color_cov(mz_cand[mask_qso], mY_cand[mask_qso],
                        mJ_cand[mask_qso], mz_err_cand[mask_qso],
                        mY_err_cand[mask_qso], mJ_err_cand[mask_qso])

# 4. Known quasars
# uki_litqso = pyfits.getdata('../resource/catalog/qso_literature_uki.fits', 1)
# redshift_litqso_uki = uki_litqso['redshift']
# fz_litqso_uki, fz_err_litqso_uki = uki_litqso['flux_z_2'], uki_litqso['flux_err_z']
# fW1_litqso_uki, fW1_err_litqso_uki = uki_litqso['flux_W1_2'], uki_litqso['flux_err_W1']
# fJ_litqso_uki, fJ_err_litqso_uki = uki_litqso['flux_J'], uki_litqso['flux_J_err']
# mz_litqso_uki, mW1_litqso_uki, mJ_litqso_uki = flux_to_mag(fz_litqso_uki), flux_to_mag(fW1_litqso_uki), flux_to_mag(fJ_litqso_uki)

# vik_litqso = pyfits.getdata('../resource/catalog/qso_literature_vik.fits', 1)
# redshift_litqso_vik = vik_litqso['redshift']
# fz_litqso_vik, fz_err_litqso_vik = vik_litqso['flux_z_2'], vik_litqso['flux_z_err']
# fW1_litqso_vik, fW1_err_litqso_vik = vik_litqso['flux_w1_2'], vik_litqso['flux_w1_err']
# fJ_litqso_vik, fJ_err_litqso_vik = vik_litqso['J_flux_aper_3p0'], vik_litqso['J_flux_aper_err_3p0']
# mz_litqso_vik, mW1_litqso_vik, mJ_litqso_vik = flux_to_mag(fz_litqso_vik), flux_to_mag(fW1_litqso_vik), flux_to_mag(fJ_litqso_uki)

# redshift_litqso = np.concatenate((redshift_litqso_uki, redshift_litqso_vik))
# mz_litqso = np.concatenate((mz_litqso_uki, mz_litqso_vik))
# mW1_litqso = np.concatenate((mW1_litqso_uki, mW1_litqso_vik))
# mJ_litqso = np.concatenate((mJ_litqso_uki, mJ_litqso_vik))
# fz_litqso = np.concatenate((fz_litqso_uki, fz_litqso_vik))
# fW1_litqso = np.concatenate((fW1_litqso_uki, fW1_litqso_vik))
# fJ_litqso = np.concatenate((fJ_litqso_uki, fJ_litqso_vik))

# 5. Quasar redshift track
sim_qso = pyfits.getdata('../resource/catalog/high_z_QSOs.fits', 1)
redshift_simqso = sim_qso['z']
redshift_bins = np.linspace(6, 8, 20)
bins_center = (redshift_bins[1:] + redshift_bins[:-1]) / 2
bins_idx = np.digitize(redshift_simqso, redshift_bins)
qso_flux = np.load('../resource/catalog/qso_flux_noisy.npy')
qso_flux_err = np.load('../resource/catalog/qso_flux_err.npy')

fJ_simqso = qso_flux[0]
fz_simqso = qso_flux[1]
fY_simqso = qso_flux[2]
fH_simqso = qso_flux[3]
fK_simqso = qso_flux[4]
fW1_simqso = qso_flux[5]
fW2_simqso = qso_flux[6]

fJ_err_simqso = qso_flux_err[0]
fz_err_simqso = qso_flux_err[1]
fY_err_simqso = qso_flux_err[2]
fH_err_simqso = qso_flux_err[3]
fK_err_simqso = qso_flux_err[4]
fW1_err_simqso = qso_flux_err[5]
fW2_err_simqso = qso_flux_err[6]

mJ_simqso = flux_to_mag(fJ_simqso, fJ_err_simqso)
mz_simqso = flux_to_mag(fz_simqso, fz_err_simqso)
mY_simqso = flux_to_mag(fY_simqso, fY_err_simqso)
mH_simqso = flux_to_mag(fH_simqso, fH_err_simqso)
mK_simqso = flux_to_mag(fK_simqso, fK_err_simqso)
mW1_simqso = flux_to_mag(fW1_simqso, fW1_err_simqso)
mW2_simqso = flux_to_mag(fW2_simqso, fW2_err_simqso)

mJ_err_simqso = fluxerr_to_magerr(fJ_err_simqso, fJ_simqso)
mz_err_simqso = fluxerr_to_magerr(fz_err_simqso, fz_simqso)
mY_err_simqso = fluxerr_to_magerr(fY_err_simqso, fY_simqso)
mH_err_simqso = fluxerr_to_magerr(fH_err_simqso, fH_simqso)
mK_err_simqso = fluxerr_to_magerr(fK_err_simqso, fK_simqso)
mW1_err_simqso = fluxerr_to_magerr(fW1_err_simqso, fW1_simqso)
mW2_err_simqso = fluxerr_to_magerr(fW2_err_simqso, fW2_simqso)

zJ_simqso = qso_flux[1] / qso_flux[0] 
YJ_simqso = qso_flux[2] / qso_flux[0]
HJ_simqso = qso_flux[3] / qso_flux[0]
KJ_simqso = qso_flux[4] / qso_flux[0]
W1J_simqso = qso_flux[5] / qso_flux[0]
W2J_simqso = qso_flux[6] / qso_flux[0]

mask_z1 = (redshift_simqso > 6.) & (redshift_simqso < 6.5) # 6.0 < z < 6.5
mask_z2 = (redshift_simqso > 6.5) & (redshift_simqso < 7) # 6.5 < z < 7.0
mask_z3 = (redshift_simqso > 7) # 7.0 < z

# 6. Dwarf spectral type track
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

syn_flux = np.load('../resource/catalog/all_syn_flux_noisy.npy')
syn_flux_err = np.load('../resource/catalog/all_syn_flux_err.npy')
syn_spectype = np.load('../resource/catalog/all_syn_spectype.npy')

fJ_syn = syn_flux[0]
fz_syn = syn_flux[1]
fY_syn = syn_flux[2]
fH_syn = syn_flux[3]
fK_syn = syn_flux[4]
fW1_syn = syn_flux[5]
fW2_syn = syn_flux[6]

fJ_err_syn = syn_flux_err[0]
fz_err_syn = syn_flux_err[1]
fY_err_syn = syn_flux_err[2]
fH_err_syn = syn_flux_err[3]
fK_err_syn = syn_flux_err[4]
fW1_err_syn = syn_flux_err[5]
fW2_err_syn = syn_flux_err[6]

mJ_syn = flux_to_mag(fJ_syn, fJ_err_syn)
mz_syn = flux_to_mag(fz_syn, fz_err_syn)
mY_syn = flux_to_mag(fY_syn, fY_err_syn)
mH_syn = flux_to_mag(fH_syn, fH_err_syn)
mK_syn = flux_to_mag(fK_syn, fK_err_syn)
mW1_syn = flux_to_mag(fW1_syn, fW1_err_syn)
mW2_syn = flux_to_mag(fW2_syn, fW2_err_syn)

mJ_err_syn = fluxerr_to_magerr(fJ_err_syn, fJ_syn)
mz_err_syn = fluxerr_to_magerr(fz_err_syn, fz_syn)
mY_err_syn = fluxerr_to_magerr(fY_err_syn, fY_syn)
mH_err_syn = fluxerr_to_magerr(fH_err_syn, fH_syn)
mK_err_syn = fluxerr_to_magerr(fK_err_syn, fK_syn)
mW1_err_syn = fluxerr_to_magerr(fW1_err_syn, fW1_syn)
mW2_err_syn = fluxerr_to_magerr(fW2_err_syn, fW2_syn)

# Summarize

all_flux = np.array([fJ_all, fz_all, fY_all, fH_all, fK_all, fW1_all, fW2_all])
all_flux_ratio = all_flux / all_flux[0, :] 
all_flux_ratio = np.delete(all_flux_ratio, 0, axis=0)
all_flux_ratio = all_flux_ratio.T
all_phot = np.array([mJ_all, mz_all, mY_all, mH_all, mK_all, mW1_all, mW2_all])
all_color = all_phot - all_phot[0, :]
all_color = np.delete(all_color, 0, axis=0)
all_color = all_color.T

cand_flux = np.array([fJ_cand, fz_cand, fY_cand, fH_cand, fK_cand, fW1_cand, fW2_cand])
cand_flux_ratio = cand_flux / cand_flux[0, :]
cand_flux_ratio = np.delete(cand_flux_ratio, 0, axis=0)
cand_flux_ratio = cand_flux_ratio.T
cand_phot = np.array([mJ_cand, mz_cand, mY_cand, mH_cand, mK_cand, mW1_cand, mW2_cand])
cand_color = cand_phot - cand_phot[0, :]
cand_color = np.delete(cand_color, 0, axis=0)
cand_color = cand_color.T

qso_flux_ratio = qso_flux / qso_flux[0, :]
qso_flux_ratio = np.delete(qso_flux_ratio, 0, axis=0)
qso_flux_ratio = qso_flux_ratio.T
qso_phot = np.array([mJ_simqso, mz_simqso, mY_simqso, mH_simqso, mK_simqso, mW1_simqso, mW2_simqso])
qso_color = qso_phot - qso_phot[0, :]
qso_color = np.delete(qso_color, 0, axis=0)
qso_color = qso_color.T

qsotrack_flux_ratio = np.zeros((len(redshift_bins)-1, 6))
qsotrack_color = np.zeros((len(redshift_bins)-1, 6))
for i in range(len(redshift_bins)-1):
    qsotrack_flux_ratio[i] = np.mean(qso_flux_ratio[bins_idx == i + 1], axis=0)
    qsotrack_color[i] = np.mean(qso_color[bins_idx == i + 1], axis=0)

# dwarf
syn_flux_ratio = syn_flux / syn_flux[0, :]
syn_flux_ratio = np.delete(syn_flux_ratio, 0, axis=0)
syn_flux_ratio = syn_flux_ratio.T
syn_phot = np.array([mJ_syn, mz_syn, mY_syn, mH_syn, mK_syn, mW1_syn, mW2_syn])
syn_color = syn_phot - syn_phot[0, :]
syn_color = np.delete(syn_color, 0, axis=0)
syn_color = syn_color.T

dwftrack_flux_ratio = np.zeros((len(spectype), 6))
dwftrack_color = np.zeros((len(spectype), 6))
for i in range(len(spectype)):
    mask = syn_spectype == spectype[i]
    dwftrack_flux_ratio[i] = np.mean(syn_flux_ratio[mask], axis=0)
    dwftrack_color[i] = np.mean(syn_color[mask], axis=0)

# Plotting

color_dwf = '#377eb8'
color_unk = '#8c510a'
color_qso = '#f781bf'
dwarf_contour_color = '#5ab4ac'
all_contour_color = 'black'

fig, ax = plt.subplots(figsize=(10, 10))

ax.scatter(cand_color[:,0][mask_Mdwarf], cand_color[:,1][mask_Mdwarf],
           s=50, marker='o', edgecolors='black', facecolors=color_dwf, 
           label=r'$\textbf{M dwarf}$', zorder=5, alpha=0.8)
ax.scatter(cand_color[:,0][mask_Ldwarf], cand_color[:,1][mask_Ldwarf],
           s=50, marker='v', edgecolors='black', facecolors=color_dwf, 
           label=r'$\textbf{L dwarf}$', zorder=5, alpha=0.8)
ax.scatter(cand_color[:,0][mask_Tdwarf], cand_color[:,1][mask_Tdwarf],
           s=50, marker='s', edgecolors='black', facecolors=color_dwf, 
           label=r'$\textbf{T dwarf}$', zorder=5, alpha=0.8)
ax.scatter(cand_color[:,0][mask_qso], cand_color[:,1][mask_qso],
           s=50, marker='o', facecolor=color_qso, edgecolor='black', 
           label=r'$\textbf{QSO}$', zorder=5)
ax.scatter(cand_color[:,0][mask_unknown], cand_color[:,1][mask_unknown],
           s=50, marker='d', facecolor=color_unk, edgecolor='black', 
           label=r'$\textbf{UNQ}$', zorder=4, alpha=0.8)              

levels = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)
ax, _ = plot_contour(all_color[:,0],
                  all_color[:,1], 
                  ax=ax, bins=50, range=((-2, 10), (-5, 10)),
                  levels=levels, 
                  smooth=.5,
                  lw=2,
                  alpha=0.4,
                  plot_datapoints=False,
                  color=all_contour_color, zorder=0)

ax, cnt_qso_z1 = plot_contour(qso_color[:,0],
                  qso_color[:,1],
                  ax=ax, bins=50, range=((-2, 10), (-5, 10)),
                  levels=levels, 
                  smooth=.5,
                  lw=2,
                  alpha=0.6,
                  plot_datapoints=False,
                  color='red', zorder=1)
# ax, cnt_qso_z1 = plot_contour(qso_color[mask_z1][:,0],
#                   qso_color[mask_z1][:,1],
#                   ax=ax, bins=50, range=((-2, 10), (-5, 10)),
#                   levels=levels, 
#                   smooth=.5,
#                   lw=2,
#                   alpha=0.6,
#                   plot_datapoints=False,
#                   color='lightsalmon', zorder=1)
# ax, cnt_qso_z2 = plot_contour(qso_color[mask_z2][:,0],
#                   qso_color[mask_z2][:,1],
#                   ax=ax, bins=50, range=((-2, 10), (-5, 10)),
#                   levels=levels, 
#                   smooth=.5,
#                   lw=2,
#                   alpha=0.6,
#                   plot_datapoints=False,
#                   color='tomato', zorder=1)
# ax, cnt_qso_z3 = plot_contour(qso_color[mask_z3][:,0],
#                   qso_color[mask_z3][:,1],
#                   ax=ax, bins=50, range=((-2, 10), (-5, 10)),
#                   levels=levels, 
#                   smooth=.5,
#                   lw=2,
#                   alpha=0.6,
#                   plot_datapoints=False,
#                   color='red', zorder=1)

ax, cnt_dwf = plot_contour(syn_color[:,0],
                  syn_color[:,1],
                  ax=ax, bins=50, range=((-2, 10), (-5, 10)),
                  levels=levels, 
                  smooth=.5,
                  lw=2,
                  alpha=0.6,
                  plot_datapoints=False,
                  color=dwarf_contour_color, zorder=2)

# qso_outer_cnt = cnt_qso_z2.collections[0].get_paths()[0]
# dwf_outer_cnt = cnt_dwf.collections[0].get_paths()[0]
# counts_qso = qso_outer_cnt.contains_points(np.vstack([syn_flux_ratio[:,0], syn_flux_ratio[:,1]]).T)
# counts_dwf = dwf_outer_cnt.contains_points(np.vstack([syn_flux_ratio[:,0], syn_flux_ratio[:,1]]).T)
# mask = np.logical_and(counts_qso, ~counts_dwf)
# ax.scatter(syn_flux_ratio[:,0][mask], syn_flux_ratio[:,1][mask], s=1)

qsotrack_x = qsotrack_color[:,0]
qsotrack_y = qsotrack_color[:,1]
mask = np.isfinite(qsotrack_x) & np.isfinite(qsotrack_y)
qsotrack_x = qsotrack_x[mask]
qsotrack_y = qsotrack_y[mask]
plot_cline(qsotrack_x, qsotrack_y, ax=ax, cmap=plt.get_cmap('autumn_r'), linewidth=2, alpha=0.5, zorder=1)
plot_cline(dwftrack_color[:,0], dwftrack_color[:,1], ax=ax, cmap=plt.get_cmap('winter_r'), linewidth=2, alpha=0.5, zorder=3)
ax.scatter(qsotrack_x, qsotrack_y, s=20, marker='o', facecolor='red', edgecolor='red', zorder=2, alpha=0.6)
ax.scatter(dwftrack_color[:,0], dwftrack_color[:,1], s=20, marker='o', facecolor=dwarf_contour_color, edgecolor=dwarf_contour_color, zorder=2, alpha=0.6)

# ax.text(-0.8, 1.5, r'${\bf z\lesssim 6.5}$', fontsize=12, color='black', zorder=5)
# ax.text(0.8, 1.5, r'${\bf z\gtrsim 6.5}$', fontsize=12, color='black', zorder=5)
# ax.text(7.2, 1.5, r'${\bf z>7}$', fontsize=12, color='black', zorder=5)

# ax.plot([-1., 0.562], [-0.261, -0.261], ls='dashed', lw=2, alpha=0.6, color='blue', zorder=2)
# ax.plot([-1., 0.562], [1.239, 1.239], ls='dashed', lw=2, alpha=0.6, color='blue', zorder=2)
# ax.plot([0.562, 0.562], [-0.261, 1.239], ls='dashed', lw=2, alpha=0.6, color='blue', zorder=2)
# ax.plot([0.562, 1.562], [1.239, 1.239], ls='dashed', lw=2, alpha=0.6, color='orange', zorder=2)
# ax.plot([0.562, 1.062], [-0.261, -0.261], ls='dashed', lw=2, alpha=0.6, color='orange', zorder=2)
# ax.plot([1.062, 1.562], [-0.261, 1.239], ls='dashed', lw=2, alpha=0.6, color='orange', zorder=2)
# ax.plot([2.062, 2.062], [-0.261, 1.239], ls='dashed', lw=2, alpha=0.6, color='red', zorder=2)
# ax.plot([2.062, 8], [-0.261, -0.261], ls='dashed', lw=2, alpha=0.6, color='red', zorder=2)
# ax.plot([2.062, 8], [1.239, 1.239], ls='dashed', lw=2, alpha=0.6, color='red', zorder=2)

ax.legend(loc='upper right', fontsize=18, frameon=True, framealpha=0.8)
ax.tick_params(axis='both', which='major', labelsize=25, width=1, size=6)
ax.tick_params(axis='both', which='minor', labelsize=25, width=1, size=3)
ax.set_xlabel(r'$z_{\rm AB}$ - $J_{\rm AB}$', fontsize=25)
ax.set_ylabel(r'$Y_{\rm AB}$ - $J_{\rm AB}$', fontsize=25)
ax.set_xlim(-1, 8)
ax.set_ylim(-1, 5.5)

ymin, ymax = ax.get_ylim()
xmin, xmax = ax.get_xlim()
xint = xmax - xmin
yint = ymax - ymin

ax.text(xmin + xint / 20, ymax - yint / 15, 
        r'$\textbf{QSOs}$', c='red',
        fontsize=20, horizontalalignment='left', weight='bold', alpha=0.8)
ax.text(xmin + xint / 20, ymax - yint / 15 - yint / 20, 
        r'$\textbf{Brown Dwarf}$ $\mathbf{(m_J=21)}$', c=dwarf_contour_color,
        fontsize=20, horizontalalignment='left', weight='bold')
ax.text(xmin + xint / 20, ymax - yint / 15 - 2 * yint / 20, 
        r'$\textbf{Observed Contaminants}$', c=all_contour_color,
        fontsize=20, horizontalalignment='left', weight='bold', alpha=0.7)

mean_cov = np.nanmean(cov_dwf_zY, axis=0)
error_cov(ax, xmax - xint / 10, ymin + yint / 10, 
          mean_cov, color=color_dwf, alpha=0.7, lw=1.5)
mean_cov = np.nanmean(cov_unk_zY, axis=0)
error_cov(ax, xmax - xint / 10, ymin + yint / 10, 
          mean_cov, color=color_unk, alpha=0.7, lw=1.5)
mean_cov = np.nanmean(cov_qso_zY, axis=0)
error_cov(ax, xmax - xint / 10, ymin + yint / 10, 
          mean_cov, color=color_qso, alpha=0.7, lw=1.5)

if args.show:
    plt.show()
else:
    # plt.savefig('cc_zY.png', dpi=200)
    plt.savefig('cc_zY.pdf', bbox_inches='tight')