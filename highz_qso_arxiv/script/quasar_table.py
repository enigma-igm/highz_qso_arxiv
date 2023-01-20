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
parser.add_argument('--meancov', action='store_true')
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

from astropy.io import fits
# cand_cat = Table(fits.getdata('../arxiv/hizqa_wildhunt.fits', 1))
cand_cat = Table(fits.getdata('../arxiv/hizqa_ric.fits', 1))
label_cand = cand_cat['label'].value
source_cand = cand_cat['source']

keys = ['name_1', 'ra_1_1', 'dec_1_1', 'redshift', 'source', \
        'flux_z', 'flux_Y', 'flux_J', 'flux_H', 'flux_K', 'flux_W1', 'flux_W2', \
        'flux_err_z', 'flux_err_Y', 'flux_err_J', 'flux_err_H', 'flux_err_K', 'flux_err_W1', 'flux_err_W2']

qso = cand_cat[(label_cand == 'QSO') & (source_cand != 'known quasar')][keys]
# convert flux to mag
qso['mag_z'] = np.round(flux_to_mag(qso['flux_z']), 2)
qso['mag_Y'] = np.round(flux_to_mag(qso['flux_Y']), 2)
qso['mag_J'] = np.round(flux_to_mag(qso['flux_J']), 2)
qso['mag_H'] = np.round(flux_to_mag(qso['flux_H']), 2)
qso['mag_K'] = np.round(flux_to_mag(qso['flux_K']), 2)
qso['mag_W1'] = np.round(flux_to_mag(qso['flux_W1']), 2)
qso['mag_W2'] = np.round(flux_to_mag(qso['flux_W2']), 2)

# TODO forget to put some error into hizqa_ric.fits

# convert flux error to mag error
# qso['mag_err_z'] = fluxerr_to_magerr(qso['flux_err_z'], qso['flux_z'])
# qso['mag_err_Y'] = fluxerr_to_magerr(qso['flux_err_Y'], qso['flux_Y'])
# qso['mag_err_J'] = fluxerr_to_magerr(qso['flux_err_J'], qso['flux_J'])
# qso['mag_err_H'] = fluxerr_to_magerr(qso['flux_err_H'], qso['flux_H'])
# qso['mag_err_K'] = fluxerr_to_magerr(qso['flux_err_K'], qso['flux_K'])
# qso['mag_err_W1'] = fluxerr_to_magerr(qso['flux_err_W1'], qso['flux_W1'])
# qso['mag_err_W2'] = fluxerr_to_magerr(qso['flux_err_W2'], qso['flux_W2'])

# qso['mag_z_with_err'] = [f'${mag:.2f} \pm {magerr:.2f}$' for mag, magerr in zip(qso['mag_z'], qso['mag_err_z'])]

# print the table into latex format
qso = qso.to_pandas()
qso = qso[['name_1', 'ra_1_1', 'dec_1_1', 'redshift', 'source', \
              'mag_z', 'mag_Y', 'mag_J', 'mag_H', 'mag_K', 'mag_W1', 'mag_W2']]
print(qso.to_latex(index=False))


