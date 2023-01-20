from astropy.io import ascii, fits
from astropy.table import Table
import pandas as pd
import numpy as np

# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('--save', action='store_true')
# parser.add_argument('--print', action='store_true')
# args = parser.parse_args()

def shorter_name(name):
    # keep only the first 4 numbers after J and +/-
    for i, nm in enumerate(name):
        if '+' in nm:
            coord = nm.split('+')
            name[i] = coord[0][:5] + '+' + coord[1][:4]
        elif '-' in nm:
            coord = nm.split('-')
            name[i] = coord[0][:5] + '-' + coord[1][:4]
    return name

def combine_jab(row, col1, col2):
    if row[col1] > 0:
        return row[col1]
    elif row[col2] > 0:
        return row[col2]
    else:
        return -1

dat = pd.read_csv('../arxiv/hizqa_catalog.csv')
wh_uki = pd.read_csv('phot/ukidssdr11las_csv/0.part')
wh_vik = pd.read_csv('phot/vikingdr5_csv/0.part')
wh_ls = pd.read_parquet('phot/ls_dr9_tractor/part.0.parquet', engine='pyarrow')

wh_uki['name'] = shorter_name(wh_uki['name'])
wh_vik['name'] = shorter_name(wh_vik['name'])
wh_ls['name'] = shorter_name(wh_ls['name'])

dat_uki = pd.merge(dat, wh_uki, on='name', how='left')
dat_vik = pd.merge(dat, wh_vik, on='name', how='left')
dat_ls = pd.merge(dat, wh_ls, on='name', how='left')

dat['flux_z'] = dat_ls['flux_z']
dat['flux_ivar_z'] = dat_ls['flux_ivar_z']
dat['flux_w1'] = dat_ls['flux_w1']
dat['flux_ivar_w1'] = dat_ls['flux_ivar_w1']
dat['flux_w2'] = dat_ls['flux_w2']
dat['flux_ivar_w2'] = dat_ls['flux_ivar_w2']

# vega to AB
# https://research.ast.cam.ac.uk/vdfs/docs/hewett-ukidss.pdf
dat['jab_uki'] = dat_uki['jAperMag3'] + 0.938
dat['j_uki_err'] = dat_uki['jAperMag3Err']
dat['yab_uki'] = dat_uki['yAperMag3'] + 0.634
dat['y_uki_err'] = dat_uki['yAperMag3Err']
dat['hab_uki'] = dat_uki['hAperMag3'] + 1.379
dat['h_uki_err'] = dat_uki['hAperMag3Err']
dat['kab_uki'] = dat_uki['kAperMag3'] + 1.900
dat['k_uki_err'] = dat_uki['kAperMag3Err']

# vega to AB: 
# http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/filter-set
dat['jab_vik'] = dat_vik['jAperMag3'] + 0.916
dat['j_vik_err'] = dat_vik['jAperMag3Err']
dat['yab_vik'] = dat_vik['yAperMag3'] + 0.600
dat['y_vik_err'] = dat_vik['yAperMag3Err']
dat['hab_vik'] = dat_vik['hAperMag3'] + 1.366
dat['h_vik_err'] = dat_vik['hAperMag3Err']
dat['kab_vik'] = dat_vik['ksAperMag3'] + 1.827
dat['k_vik_err'] = dat_vik['ksAperMag3Err']

dat['jab'] = dat.apply(combine_jab, axis=1, col1='jab_uki', col2='jab_vik')
dat['j_err'] = dat.apply(combine_jab, axis=1, col1='j_uki_err', col2='j_vik_err')
dat['yab'] = dat.apply(combine_jab, axis=1, col1='yab_uki', col2='yab_vik')
dat['y_err'] = dat.apply(combine_jab, axis=1, col1='y_uki_err', col2='y_vik_err')
dat['hab'] = dat.apply(combine_jab, axis=1, col1='hab_uki', col2='hab_vik')
dat['h_err'] = dat.apply(combine_jab, axis=1, col1='h_uki_err', col2='h_vik_err')
dat['kab'] = dat.apply(combine_jab, axis=1, col1='kab_uki', col2='kab_vik')
dat['k_err'] = dat.apply(combine_jab, axis=1, col1='k_uki_err', col2='k_vik_err')
dat['JAB'] = np.round(dat['jab'], 2)
dat['YAB'] = np.round(dat['yab'], 2)
dat['HAB'] = np.round(dat['hab'], 2)
dat['KAB'] = np.round(dat['kab'], 2)
dat = dat.drop(columns=['jab_uki', 'jab_vik', 'jab'])
dat = dat.drop(columns=['j_uki_err', 'j_vik_err'])
dat = dat.drop(columns=['yab_uki', 'yab_vik', 'yab'])
dat = dat.drop(columns=['y_uki_err', 'y_vik_err'])
dat = dat.drop(columns=['hab_uki', 'hab_vik', 'hab'])
dat = dat.drop(columns=['h_uki_err', 'h_vik_err'])
dat = dat.drop(columns=['kab_uki', 'kab_vik', 'kab'])
dat = dat.drop(columns=['k_uki_err', 'k_vik_err'])

dat.to_csv('../arxiv/hizqa_catalog_wildhunt.csv', index=False)
dat = Table.from_pandas(dat)
dat.write('../arxiv/hizqa_wildhunt.fits', overwrite=True)

# if args.save:
#     dat.write('hizqa_catalog.csv')
# elif args.print:
#     latex = dat['name', 'ra', 'dec', 'JAB', 'label', 'spectral type', 'redshift', 'source', 'run'].to_pandas().to_latex(index=False)
#     print(latex)