from math import comb
from astropy.io import ascii, fits
from astropy.table import Table
import pandas as pd
import numpy as np

def combine_jab(row, col1, col2):
    if row[col1] > 0:
        return row[col1]
    elif row[col2] > 0:
        return row[col2]
    else:
        return -1

dat = Table(fits.getdata('../arxiv/hizqa_ric.fits', 1))
dat = dat.to_pandas()

dat['flux_z'] = dat.apply(combine_jab, axis=1, col1='flux_z_1', col2='flux_z')
dat['flux_err_z'] = dat.apply(combine_jab, axis=1, col1='flux_z_err', col2='flux_err_z')
dat['flux_J'] = dat.apply(combine_jab, axis=1, col1='flux_J', col2='J_flux_aper_3p0')
dat['flux_err_J'] = dat.apply(combine_jab, axis=1, col1='flux_J_err', col2='J_flux_aper_err_3p0')
dat['flux_H'] = dat.apply(combine_jab, axis=1, col1='H_flux_aper_3p0_1', col2='H_flux_aper_3p0_2')
dat['flux_err_H'] = dat.apply(combine_jab, axis=1, col1='H_flux_aper_err_3p0_1', col2='H_flux_aper_err_3p0_2')
dat['flux_K'] = dat.apply(combine_jab, axis=1, col1='K_flux_aper_3p0_1', col2='K_flux_aper_3p0_2')
dat['flux_err_K'] = dat.apply(combine_jab, axis=1, col1='K_flux_aper_err_3p0_1', col2='K_flux_aper_err_3p0_2')
dat['flux_Y'] = dat.apply(combine_jab, axis=1, col1='Y_flux_aper_3p0_1', col2='Y_flux_aper_3p0_2')
dat['flux_err_Y'] = dat.apply(combine_jab, axis=1, col1='Y_flux_aper_err_3p0_1', col2='Y_flux_aper_err_3p0_2')
dat['flux_W1'] = dat.apply(combine_jab, axis=1, col1='flux_w1_2', col2='flux_W1')
dat['flux_err_W1'] = dat.apply(combine_jab, axis=1, col1='flux_w1_err', col2='flux_err_W1')
dat['flux_W2'] = dat.apply(combine_jab, axis=1, col1='flux_w2_2', col2='flux_W2')
dat['flux_err_W2'] = dat.apply(combine_jab, axis=1, col1='flux_w2_err', col2='flux_err_W2')

dat = Table.from_pandas(dat)
dat.write('../arxiv/hizqa_ric.fits', overwrite=True)