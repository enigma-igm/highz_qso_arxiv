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

def flux_to_mag(flux):
    mask = flux > 0
    mag = np.zeros_like(flux)
    mag[mask] = -2.5 * np.log10(flux[mask]) + 22.5
    mag[~mask] = -1
    return mag
def mag_to_flux(mag):
    mask = mag > 0 # dangerous?!
    flux = np.zeros_like(mag)
    flux[mask] = 10**(-0.4 * (mag[mask] - 22.5))
    flux[~mask] = -1
    return flux
def magerr_to_fluxerr(magerr, flux):
    mask = magerr > 0
    fluxerr = np.zeros_like(magerr)
    fluxerr[mask] = flux[mask] * np.log(10) / 2.5 * magerr[mask]
    fluxerr[~mask] = -1
    return fluxerr
def fluxerr_to_magerr(fluxerr, flux):
    mask = fluxerr > 0
    magerr = np.zeros_like(fluxerr)
    magerr[mask] = 2.5 / np.log(10) * fluxerr[mask] / flux[mask]
    magerr[~mask] = -1
    return magerr
def flux_ratio(color):
    return 10**(-color/2.5)
def flux_ratio_cov(X, Y, Z, Xerr, Yerr, Zerr):
    # X/Z vs. Y/Z
    cov = np.zeros((len(X), 2, 2))
    cov[:, 0, 0] = (X / Z)**2 * ((Xerr / X)**2 + (Zerr / Z)**2)
    cov[:, 1, 1] = (Y / Z)**2 * ((Yerr / Y)**2 + (Zerr / Z)**2)
    cov[:, 0, 1] = cov[:, 1, 0] = X * Y * Zerr**2 / Z**4
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

# TODO: add non-detection
mz_uki, mW1_uki, mJ_uki = flux_to_mag(fz_uki), flux_to_mag(fW1_uki), flux_to_mag(fJ_uki)
mW2_uki, mY_uki, mH_uki, mK_uki = flux_to_mag(fW2_uki), flux_to_mag(fY_uki), flux_to_mag(fH_uki), flux_to_mag(fK_uki)

# all viking targets
vik_all = pyfits.getdata('../resource/catalog/VIKING_catalog_clean_nobright.fits', 1)
fz_vik, fz_err_vik = vik_all['flux_z'], vik_all['flux_z_err']
fW1_vik, fW1_err_vik = vik_all['flux_w1'], vik_all['flux_w1_err']
fW2_vik, fW2_err_vik = vik_all['flux_w2'], vik_all['flux_w2_err']
fY_vik, fY_err_vik = vik_all['Y_flux_aper_3p0'], vik_all['Y_flux_aper_err_3p0']
fJ_vik, fJ_err_vik = vik_all['J_flux_aper_3p0'], vik_all['J_flux_aper_err_3p0']
fH_vik, fH_err_vik = vik_all['H_flux_aper_3p0'], vik_all['H_flux_aper_err_3p0']
fK_vik, fK_err_vik = vik_all['K_flux_aper_3p0'], vik_all['K_flux_aper_err_3p0']


# TODO: add non-detection
mz_vik, mW1_vik, mJ_vik = flux_to_mag(fz_vik), flux_to_mag(fW1_vik), flux_to_mag(fJ_vik)
mW2_vik, mY_vik, mH_vik, mK_vik = flux_to_mag(fW2_vik), flux_to_mag(fY_vik), flux_to_mag(fH_vik), flux_to_mag(fK_vik)

mz_all = np.concatenate([mz_uki, mz_vik])
mW1_all = np.concatenate([mW1_uki, mW1_vik])
mW2_all = np.concatenate([mW2_uki, mW2_vik])
mJ_all = np.concatenate([mJ_uki, mJ_vik])
mask = np.isfinite(mz_all) & np.isfinite(mW1_all) & np.isfinite(mJ_all)
fz_all = np.concatenate([fz_uki, fz_vik])
fW1_all = np.concatenate([fW1_uki, fW1_vik])
fW2_all = np.concatenate([fW2_uki, fW2_vik])
fJ_all = np.concatenate([fJ_uki, fJ_vik])
fY_all = np.concatenate([fY_uki, fY_vik])
fH_all = np.concatenate([fH_uki, fH_vik])
fK_all = np.concatenate([fK_uki, fK_vik])

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

cov_zY = flux_ratio_cov(fz_cand[mask_detected], fY_cand[mask_detected], fJ_cand[mask_detected],
                        fz_err_cand[mask_detected], fY_err_cand[mask_detected], fJ_err_cand[mask_detected])
cov_dwf_zY = flux_ratio_cov(fz_cand[mask_dwarf], fY_cand[mask_dwarf], fJ_cand[mask_dwarf],
                            fz_err_cand[mask_dwarf], fY_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_zY = flux_ratio_cov(fz_cand[mask_unknown], fY_cand[mask_unknown], fJ_cand[mask_unknown],
                            fz_err_cand[mask_unknown], fY_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_zY = flux_ratio_cov(fz_cand[mask_qso], fY_cand[mask_qso], fJ_cand[mask_qso],
                            fz_err_cand[mask_qso], fY_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_zH = flux_ratio_cov(fz_cand[mask_detected], fH_cand[mask_detected], fJ_cand[mask_detected],
                        fz_err_cand[mask_detected], fH_err_cand[mask_detected], fJ_err_cand[mask_detected])
cov_dwf_zH = flux_ratio_cov(fz_cand[mask_dwarf], fH_cand[mask_dwarf], fJ_cand[mask_dwarf],
                            fz_err_cand[mask_dwarf], fH_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_zH = flux_ratio_cov(fz_cand[mask_unknown], fH_cand[mask_unknown], fJ_cand[mask_unknown],
                            fz_err_cand[mask_unknown], fH_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_zH = flux_ratio_cov(fz_cand[mask_qso], fH_cand[mask_qso], fJ_cand[mask_qso],
                            fz_err_cand[mask_qso], fH_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_zK = flux_ratio_cov(fz_cand[mask_detected], fK_cand[mask_detected], fJ_cand[mask_detected],
                        fz_err_cand[mask_detected], fK_err_cand[mask_detected], fJ_err_cand[mask_detected])
cov_dwf_zK = flux_ratio_cov(fz_cand[mask_dwarf], fK_cand[mask_dwarf], fJ_cand[mask_dwarf],
                            fz_err_cand[mask_dwarf], fK_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_zK = flux_ratio_cov(fz_cand[mask_unknown], fK_cand[mask_unknown], fJ_cand[mask_unknown],
                            fz_err_cand[mask_unknown], fK_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_zK = flux_ratio_cov(fz_cand[mask_qso], fK_cand[mask_qso], fJ_cand[mask_qso],
                            fz_err_cand[mask_qso], fK_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_zW1 = flux_ratio_cov(fz_cand[mask_detected], fW1_cand[mask_detected], fJ_cand[mask_detected],
                         fz_err_cand[mask_detected], fW1_err_cand[mask_detected], fJ_err_cand[mask_detected])
cov_dwf_zW1 = flux_ratio_cov(fz_cand[mask_dwarf], fW1_cand[mask_dwarf], fJ_cand[mask_dwarf],
                             fz_err_cand[mask_dwarf], fW1_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_zW1 = flux_ratio_cov(fz_cand[mask_unknown], fW1_cand[mask_unknown], fJ_cand[mask_unknown],
                             fz_err_cand[mask_unknown], fW1_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_zW1 = flux_ratio_cov(fz_cand[mask_qso], fW1_cand[mask_qso], fJ_cand[mask_qso],
                             fz_err_cand[mask_qso], fW1_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_zW2 = flux_ratio_cov(fz_cand[mask_detected], fW2_cand[mask_detected], fJ_cand[mask_detected],
                         fz_err_cand[mask_detected], fW2_err_cand[mask_detected], fJ_err_cand[mask_detected])
cov_dwf_zW2 = flux_ratio_cov(fz_cand[mask_dwarf], fW2_cand[mask_dwarf], fJ_cand[mask_dwarf],
                             fz_err_cand[mask_dwarf], fW2_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_zW2 = flux_ratio_cov(fz_cand[mask_unknown], fW2_cand[mask_unknown], fJ_cand[mask_unknown],
                             fz_err_cand[mask_unknown], fW2_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_zW2 = flux_ratio_cov(fz_cand[mask_qso], fW2_cand[mask_qso], fJ_cand[mask_qso],
                             fz_err_cand[mask_qso], fW2_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_YH = flux_ratio_cov(fY_cand[mask_detected], fH_cand[mask_detected], fJ_cand[mask_detected],
                        fY_err_cand[mask_detected], fH_err_cand[mask_detected], fJ_err_cand[mask_detected])
cov_dwf_YH = flux_ratio_cov(fY_cand[mask_dwarf], fH_cand[mask_dwarf], fJ_cand[mask_dwarf],
                            fY_err_cand[mask_dwarf], fH_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_YH = flux_ratio_cov(fY_cand[mask_unknown], fH_cand[mask_unknown], fJ_cand[mask_unknown],
                            fY_err_cand[mask_unknown], fH_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_YH = flux_ratio_cov(fY_cand[mask_qso], fH_cand[mask_qso], fJ_cand[mask_qso],
                            fY_err_cand[mask_qso], fH_err_cand[mask_qso], fJ_err_cand[mask_qso])
                            
cov_YK = flux_ratio_cov(fY_cand[mask_detected], fK_cand[mask_detected], fJ_cand[mask_detected],
                        fY_err_cand[mask_detected], fK_err_cand[mask_detected], fJ_err_cand[mask_detected]) 
cov_dwf_YK = flux_ratio_cov(fY_cand[mask_dwarf], fK_cand[mask_dwarf], fJ_cand[mask_dwarf],
                            fY_err_cand[mask_dwarf], fK_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])  
cov_unk_YK = flux_ratio_cov(fY_cand[mask_unknown], fK_cand[mask_unknown], fJ_cand[mask_unknown],
                            fY_err_cand[mask_unknown], fK_err_cand[mask_unknown], fJ_err_cand[mask_unknown])  
cov_qso_YK = flux_ratio_cov(fY_cand[mask_qso], fK_cand[mask_qso], fJ_cand[mask_qso],
                            fY_err_cand[mask_qso], fK_err_cand[mask_qso], fJ_err_cand[mask_qso])    

cov_YW1 = flux_ratio_cov(fY_cand[mask_detected], fW1_cand[mask_detected], fJ_cand[mask_detected],
                         fY_err_cand[mask_detected], fW1_err_cand[mask_detected], fJ_err_cand[mask_detected])
cov_dwf_YW1 = flux_ratio_cov(fY_cand[mask_dwarf], fW1_cand[mask_dwarf], fJ_cand[mask_dwarf],
                             fY_err_cand[mask_dwarf], fW1_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_YW1 = flux_ratio_cov(fY_cand[mask_unknown], fW1_cand[mask_unknown], fJ_cand[mask_unknown],
                             fY_err_cand[mask_unknown], fW1_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_YW1 = flux_ratio_cov(fY_cand[mask_qso], fW1_cand[mask_qso], fJ_cand[mask_qso],
                             fY_err_cand[mask_qso], fW1_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_YW2 = flux_ratio_cov(fY_cand[mask_detected], fW2_cand[mask_detected], fJ_cand[mask_detected],
                         fY_err_cand[mask_detected], fW2_err_cand[mask_detected], fJ_err_cand[mask_detected])    
cov_dwf_YW2 = flux_ratio_cov(fY_cand[mask_dwarf], fW2_cand[mask_dwarf], fJ_cand[mask_dwarf],
                             fY_err_cand[mask_dwarf], fW2_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_YW2 = flux_ratio_cov(fY_cand[mask_unknown], fW2_cand[mask_unknown], fJ_cand[mask_unknown],
                             fY_err_cand[mask_unknown], fW2_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_YW2 = flux_ratio_cov(fY_cand[mask_qso], fW2_cand[mask_qso], fJ_cand[mask_qso],
                             fY_err_cand[mask_qso], fW2_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_HK = flux_ratio_cov(fH_cand[mask_detected], fK_cand[mask_detected], fJ_cand[mask_detected],
                        fH_err_cand[mask_detected], fK_err_cand[mask_detected], fJ_err_cand[mask_detected])       
cov_dwf_HK = flux_ratio_cov(fH_cand[mask_dwarf], fK_cand[mask_dwarf], fJ_cand[mask_dwarf],
                            fH_err_cand[mask_dwarf], fK_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_HK = flux_ratio_cov(fH_cand[mask_unknown], fK_cand[mask_unknown], fJ_cand[mask_unknown],
                            fH_err_cand[mask_unknown], fK_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_HK = flux_ratio_cov(fH_cand[mask_qso], fK_cand[mask_qso], fJ_cand[mask_qso],
                            fH_err_cand[mask_qso], fK_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_HW1 = flux_ratio_cov(fH_cand[mask_detected], fW1_cand[mask_detected], fJ_cand[mask_detected],
                         fH_err_cand[mask_detected], fW1_err_cand[mask_detected], fJ_err_cand[mask_detected])        
cov_dwf_HW1 = flux_ratio_cov(fH_cand[mask_dwarf], fW1_cand[mask_dwarf], fJ_cand[mask_dwarf],
                             fH_err_cand[mask_dwarf], fW1_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_HW1 = flux_ratio_cov(fH_cand[mask_unknown], fW1_cand[mask_unknown], fJ_cand[mask_unknown],
                             fH_err_cand[mask_unknown], fW1_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_HW1 = flux_ratio_cov(fH_cand[mask_qso], fW1_cand[mask_qso], fJ_cand[mask_qso],
                             fH_err_cand[mask_qso], fW1_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_HW2 = flux_ratio_cov(fH_cand[mask_detected], fW2_cand[mask_detected], fJ_cand[mask_detected],
                         fH_err_cand[mask_detected], fW2_err_cand[mask_detected], fJ_err_cand[mask_detected])
cov_dwf_HW2 = flux_ratio_cov(fH_cand[mask_dwarf], fW2_cand[mask_dwarf], fJ_cand[mask_dwarf],
                             fH_err_cand[mask_dwarf], fW2_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_HW2 = flux_ratio_cov(fH_cand[mask_unknown], fW2_cand[mask_unknown], fJ_cand[mask_unknown],
                             fH_err_cand[mask_unknown], fW2_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_HW2 = flux_ratio_cov(fH_cand[mask_qso], fW2_cand[mask_qso], fJ_cand[mask_qso],
                             fH_err_cand[mask_qso], fW2_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_KW1 = flux_ratio_cov(fK_cand[mask_detected], fW1_cand[mask_detected], fJ_cand[mask_detected],
                         fK_err_cand[mask_detected], fW1_err_cand[mask_detected], fJ_err_cand[mask_detected])    
cov_dwf_KW1 = flux_ratio_cov(fK_cand[mask_dwarf], fW1_cand[mask_dwarf], fJ_cand[mask_dwarf],
                             fK_err_cand[mask_dwarf], fW1_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_KW1 = flux_ratio_cov(fK_cand[mask_unknown], fW1_cand[mask_unknown], fJ_cand[mask_unknown],
                             fK_err_cand[mask_unknown], fW1_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_KW1 = flux_ratio_cov(fK_cand[mask_qso], fW1_cand[mask_qso], fJ_cand[mask_qso],
                             fK_err_cand[mask_qso], fW1_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_KW2 = flux_ratio_cov(fK_cand[mask_detected], fW2_cand[mask_detected], fJ_cand[mask_detected],
                         fK_err_cand[mask_detected], fW2_err_cand[mask_detected], fJ_err_cand[mask_detected])
cov_dwf_KW2 = flux_ratio_cov(fK_cand[mask_dwarf], fW2_cand[mask_dwarf], fJ_cand[mask_dwarf],
                             fK_err_cand[mask_dwarf], fW2_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_KW2 = flux_ratio_cov(fK_cand[mask_unknown], fW2_cand[mask_unknown], fJ_cand[mask_unknown],
                             fK_err_cand[mask_unknown], fW2_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_KW2 = flux_ratio_cov(fK_cand[mask_qso], fW2_cand[mask_qso], fJ_cand[mask_qso],
                             fK_err_cand[mask_qso], fW2_err_cand[mask_qso], fJ_err_cand[mask_qso])

cov_W1W2 = flux_ratio_cov(fW1_cand[mask_detected], fW2_cand[mask_detected], fJ_cand[mask_detected],
                          fW1_err_cand[mask_detected], fW2_err_cand[mask_detected], fJ_err_cand[mask_detected])        
cov_dwf_W1W2 = flux_ratio_cov(fW1_cand[mask_dwarf], fW2_cand[mask_dwarf], fJ_cand[mask_dwarf],
                              fW1_err_cand[mask_dwarf], fW2_err_cand[mask_dwarf], fJ_err_cand[mask_dwarf])
cov_unk_W1W2 = flux_ratio_cov(fW1_cand[mask_unknown], fW2_cand[mask_unknown], fJ_cand[mask_unknown],
                              fW1_err_cand[mask_unknown], fW2_err_cand[mask_unknown], fJ_err_cand[mask_unknown])
cov_qso_W1W2 = flux_ratio_cov(fW1_cand[mask_qso], fW2_cand[mask_qso], fJ_cand[mask_qso],
                              fW1_err_cand[mask_qso], fW2_err_cand[mask_qso], fJ_err_cand[mask_qso])

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
redshift_bins = np.linspace(6, 8, 50)
bins_center = (redshift_bins[1:] + redshift_bins[:-1]) / 2
bins_idx = np.digitize(redshift_simqso, redshift_bins)
qso_flux = np.load('../resource/catalog/qso_flux_noisy.npy')
# qso_flux = np.array([fJ_simqso, fz_simqso, fY_simqso, fH_simqso, fK_simqso, fW1_simqso, fW2_simqso])

zJ_simqso = qso_flux[1] / qso_flux[0] 
YJ_simqso = qso_flux[2] / qso_flux[0]
HJ_simqso = qso_flux[3] / qso_flux[0]
KJ_simqso = qso_flux[4] / qso_flux[0]
W1J_simqso = qso_flux[5] / qso_flux[0]
W2J_simqso = qso_flux[6] / qso_flux[0]

zJ_simqso_mean = np.zeros(len(redshift_bins) - 1)
YJ_simqso_mean = np.zeros(len(redshift_bins) - 1)
HJ_simqso_mean = np.zeros(len(redshift_bins) - 1)
KJ_simqso_mean = np.zeros(len(redshift_bins) - 1)
W1J_simqso_mean = np.zeros(len(redshift_bins) - 1)
W2J_simqso_mean = np.zeros(len(redshift_bins) - 1)
for i in range(len(redshift_bins) - 1):
    zJ_simqso_mean[i] = np.median(zJ_simqso[bins_idx == i + 1])
    YJ_simqso_mean[i] = np.median(YJ_simqso[bins_idx == i + 1])
    HJ_simqso_mean[i] = np.median(HJ_simqso[bins_idx == i + 1])
    KJ_simqso_mean[i] = np.median(KJ_simqso[bins_idx == i + 1])
    W1J_simqso_mean[i] = np.median(W1J_simqso[bins_idx == i + 1])
    W2J_simqso_mean[i] = np.median(W2J_simqso[bins_idx == i + 1])

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

# Summarize

all_flux = np.array([fJ_all, fz_all, fY_all, fH_all, fK_all, fW1_all, fW2_all])
all_flux_ratio = all_flux / all_flux[0, :] 
all_flux_ratio = np.delete(all_flux_ratio, 0, axis=0)
all_flux_ratio = all_flux_ratio.T

cand_flux = np.array([fJ_cand, fz_cand, fY_cand, fH_cand, fK_cand, fW1_cand, fW2_cand])
cand_flux_ratio = cand_flux / cand_flux[0, :]
cand_flux_ratio = np.delete(cand_flux_ratio, 0, axis=0)
cand_flux_ratio = cand_flux_ratio.T

dwftrack_flux_ratio = np.array([fzJ, fYJ, fHJ, fKJ, fW1J, fW2J])
dwftrack_flux_ratio = dwftrack_flux_ratio.T

qso_flux_ratio = qso_flux / qso_flux[0, :]
qso_flux_ratio = np.delete(qso_flux_ratio, 0, axis=0)
qso_flux_ratio = qso_flux_ratio.T

qsotrack_flux = np.array([zJ_simqso_mean, YJ_simqso_mean, HJ_simqso_mean, KJ_simqso_mean, W1J_simqso_mean, W2J_simqso_mean])
qsotrack_flux = qsotrack_flux.T

syn_flux = np.load('../resource/catalog/all_syn_flux_noisy.npy')
syn_flux_ratio = syn_flux / syn_flux[0, :]
syn_flux_ratio = np.delete(syn_flux_ratio, 0, axis=0)
syn_flux_ratio = syn_flux_ratio.T

# Plotting

dwarf_color = '#377eb8'
unknown_color = '#8c510a'
qso_color = '#f781bf'
dwarf_contour_color = '#5ab4ac'
all_contour_color = 'black'

fig, ax = plt.subplots(figsize=(10, 10))

# ax.text(0.5, -0.6, r'$\mathbf{6<z<6.5}$ $\textbf{QSOs}$', c='lightsalmon',
#               fontsize=20, horizontalalignment='left', weight='bold')
# ax.text(0.5, -0.7, r'$\mathbf{6.5<z<7}$ $\textbf{QSOs}$', c='tomato',
#               fontsize=20, horizontalalignment='left', weight='bold')
# ax.text(0.5, -0.8, r'$\mathbf{z>7}$ $\textbf{QSOs}$', c='red',
#               fontsize=20, horizontalalignment='left', weight='bold')
# ax.text(0.5, -0.9, r'$\textbf{Brown Dwarf}$', c=dwarf_contour_color,
#                 fontsize=20, horizontalalignment='left', weight='bold')

ax.text(-0.04, -0.7, r'$\textbf{QSOs}$', c='red',
              fontsize=20, horizontalalignment='left', weight='bold', alpha=0.8)
ax.text(-0.04, -0.8, r'$\textbf{Brown Dwarf}$ ($\mathbf{m_J=21}$)', c=dwarf_contour_color,
                fontsize=20, horizontalalignment='left', weight='bold')
ax.text(-0.04, -0.9, r'$\textbf{Observed Contaminants}$', c=all_contour_color,
              fontsize=20, horizontalalignment='left', weight='bold', alpha=0.7)

ax.scatter(cand_flux[1][mask_Mdwarf] / cand_flux[0][mask_Mdwarf], cand_flux[2][mask_Mdwarf] / cand_flux[0][mask_Mdwarf], 
            s=50, marker='o', edgecolors='black', facecolors=dwarf_color, label=r'$\textbf{M dwarf}$', zorder=5, alpha=0.8)
ax.scatter(cand_flux[1][mask_Ldwarf] / cand_flux[0][mask_Ldwarf], cand_flux[2][mask_Ldwarf] / cand_flux[0][mask_Ldwarf], 
            s=50, marker='v', edgecolors='black', facecolors=dwarf_color, label=r'$\textbf{L dwarf}$', zorder=5, alpha=0.8)
ax.scatter(cand_flux[1][mask_Tdwarf] / cand_flux[0][mask_Tdwarf], cand_flux[2][mask_Tdwarf] / cand_flux[0][mask_Tdwarf], 
            s=50, marker='s', edgecolors='black', facecolors=dwarf_color, label=r'$\textbf{T dwarf}$', zorder=5, alpha=0.8)
ax.scatter(cand_flux[1][mask_qso] / cand_flux[0][mask_qso], cand_flux[2][mask_qso] / cand_flux[0][mask_qso],
            s=50, marker='o', facecolor=qso_color, edgecolor='black', label=r'$\textbf{QSO}$', zorder=5)
ax.scatter(cand_flux[1][mask_unknown] / cand_flux[0][mask_unknown], cand_flux[2][mask_unknown] / cand_flux[0][mask_unknown],
            s=50, marker='d', facecolor=unknown_color, edgecolor='black', label=r'$\textbf{UNQ}$', zorder=4, alpha=0.8)

mean_cov = np.nanmean(cov_dwf_zY, axis=0)
error_cov(ax, 1.25, -0.7, mean_cov, color=dwarf_color, alpha=0.7, lw=1.5)
mean_cov = np.nanmean(cov_unk_zY, axis=0)
error_cov(ax, 1.25, -0.7, mean_cov, color=unknown_color, alpha=0.7, lw=1.5)
mean_cov = np.nanmean(cov_qso_zY, axis=0)
error_cov(ax, 1.25, -0.7, mean_cov, color=qso_color, alpha=0.7, lw=1.5)

levels = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)
ax, _ = plot_contour(all_flux_ratio[:,0],
                  all_flux_ratio[:,1], 
                  ax=ax, bins=50, range=((-0.5,1.), (-0.7,1.4)),
                  levels=levels, 
                  smooth=.5,
                  lw=2,
                  alpha=0.4,
                  plot_datapoints=False,
                  color=all_contour_color, zorder=0)

ax, cnt_qso_z1 = plot_contour(qso_flux_ratio[:,0],
                  qso_flux_ratio[:,1],
                  ax=ax, bins=50, range=((-0.5,1.), (-0.7,1.4)),
                  levels=levels, 
                  smooth=.5,
                  lw=2,
                  alpha=0.6,
                  plot_datapoints=False,
                  color='red', zorder=1)
# ax, cnt_qso_z1 = plot_contour(qso_flux_ratio[mask_z1][:,0],
#                   qso_flux_ratio[mask_z1][:,1],
#                   ax=ax, bins=50, range=((-0.5,1.), (-0.7,1.4)),
#                   levels=levels, 
#                   smooth=.5,
#                   lw=2,
#                   alpha=0.6,
#                   plot_datapoints=False,
#                   color='lightsalmon', zorder=1)
# ax, cnt_qso_z2 = plot_contour(qso_flux_ratio[mask_z2][:,0],
#                   qso_flux_ratio[mask_z2][:,1],
#                   ax=ax, bins=50, range=((-0.5,1.), (-0.7,1.4)),
#                   levels=levels, 
#                   smooth=.5,
#                   lw=2,
#                   alpha=0.6,
#                   plot_datapoints=False,
#                   color='tomato', zorder=1)
# ax, cnt_qso_z3 = plot_contour(qso_flux_ratio[mask_z3][:,0],
#                   qso_flux_ratio[mask_z3][:,1],
#                   ax=ax, bins=50, range=((-0.5,1.), (-0.7,1.4)),
#                   levels=levels, 
#                   smooth=.5,
#                   lw=2,
#                   alpha=0.6,
#                   plot_datapoints=False,
#                   color='red', zorder=1)

ax, cnt_dwf = plot_contour(syn_flux_ratio[:,0],
                  syn_flux_ratio[:,1],
                  ax=ax, bins=50, range=((-0.5,1.), (-0.7,1.4)),
                  levels=levels, 
                  smooth=.5,
                  lw=2,
                  alpha=0.6,
                  plot_datapoints=False,
                  color=dwarf_contour_color, zorder=2)

plot_cline(qsotrack_flux[:,0], qsotrack_flux[:,1], ax=ax, cmap=plt.get_cmap('autumn_r'), linewidth=2, alpha=0.5, zorder=1)
plot_cline(dwftrack_flux_ratio[:,0], dwftrack_flux_ratio[:,1], ax=ax, cmap=plt.get_cmap('winter_r'), linewidth=2, alpha=0.5, zorder=3)
ax.scatter(qsotrack_flux[:,0], qsotrack_flux[:,1], s=20, marker='o', facecolor='red', edgecolor='red', zorder=2, alpha=0.6)
ax.scatter(dwftrack_flux_ratio[:,0], dwftrack_flux_ratio[:,1], s=20, marker='o', facecolor=dwarf_contour_color, edgecolor=dwarf_contour_color, zorder=2, alpha=0.6)

# ax.plot(flux_ratio(np.array([-1., 0.562])), 1/flux_ratio(np.array([-0.261, -0.261])), ls='dashed', lw=2, alpha=0.6, color='blue', zorder=2)
# ax.plot(flux_ratio(np.array([-1., 0.562])), 1/flux_ratio(np.array([1.239, 1.239])), ls='dashed', lw=2, alpha=0.6, color='blue', zorder=2)
# ax.plot(flux_ratio(np.array([0.562, 0.562])), 1/flux_ratio(np.array([-0.261, 1.239])), ls='dashed', lw=2, alpha=0.6, color='blue', zorder=2)
# ax.plot(flux_ratio(np.array([0.562, 1.562])), 1/flux_ratio(np.array([1.239, 1.239])), ls='dashed', lw=2, alpha=0.6, color='orange', zorder=2)
# ax.plot(flux_ratio(np.array([0.562, 1.062])), 1/flux_ratio(np.array([-0.261, -0.261])), ls='dashed', lw=2, alpha=0.6, color='orange', zorder=2)
# ax.plot(flux_ratio(np.array([1.062, 1.562])), 1/flux_ratio(np.array([-0.261, 1.239])), ls='dashed', lw=2, alpha=0.6, color='orange', zorder=2)
# ax.plot(flux_ratio(np.array([2.062, 2.062])), 1/flux_ratio(np.array([-0.261, 1.239])), ls='dashed', lw=2, alpha=0.6, color='red', zorder=2)
# x = flux_ratio(np.array([2.062, 8]))
# x[1] = -0.1 
# ax.plot(x, 1/flux_ratio(np.array([-0.261, -0.261])), ls='dashed', lw=2, alpha=0.6, color='red', zorder=2)
# ax.plot(x, 1/flux_ratio(np.array([1.239, 1.239])), ls='dashed', lw=2, alpha=0.6, color='red', zorder=2)
# ax.text(1.3, 2.9, r'${\bf z\lesssim 6.5}$', fontsize=12, color='black', zorder=5)
# ax.text(0.42, 2.9, r'${\bf z\gtrsim 6.5}$', fontsize=12, color='black', zorder=5)
# ax.text(-0.02, 2.9, r'${\bf z>7}$', fontsize=12, color='black', zorder=5)

ax.legend(loc='upper right', fontsize=18, frameon=True, framealpha=0.8)
ax.tick_params(axis='both', which='major', labelsize=25, width=1, size=6)
ax.tick_params(axis='both', which='minor', labelsize=25, width=1, size=3)
ax.set_xlabel(r'$f_{z} / f_{J}$', fontsize=25)
ax.set_ylabel(r'$f_{Y} / f_{J}$', fontsize=25)

ax.set_xlim(-0.1, 1.5)
ax.set_ylim(-1,1.6)

if args.show:
    plt.show()
else:
    # plt.savefig('flux_ratio_zY.png', dpi=200)
    plt.savefig('flux_ratio_zY.pdf', bbox_inches='tight')