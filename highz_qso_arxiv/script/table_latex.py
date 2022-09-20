from astropy.io import ascii, fits
from astropy.table import Table
import numpy as np

dat = Table(fits.getdata('../arxiv/hizqa_catalog.fits', 1))

name = dat['name']
# keep only the first 4 numbers after J and +/-
for i, nm in enumerate(name):
    if '+' in nm:
        coord = nm.split('+')
        name[i] = coord[0][:5] + '+' + coord[1][:4]
    elif '-' in nm:
        coord = nm.split('-')
        name[i] = coord[0][:5] + '-' + coord[1][:4]
dat['name'] = name
dat['ra'] = np.round(dat['ra'], 4)
dat['dec'] = np.round(dat['dec'], 4)
dat['JAB'] = np.round(dat['JAB'], 2)

latex = dat['name', 'ra', 'dec', 'JAB', 'label', 'source'].to_pandas().to_latex(index=False)
print(latex)
