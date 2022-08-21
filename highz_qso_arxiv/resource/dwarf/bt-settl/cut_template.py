import os
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

from astropy.io import ascii
from astropy.table import Table
from scipy.interpolate import interp1d
from highz_qso_arxiv.util import get_project_root
from highz_qso_arxiv.resource.filters import ukirt_J, wise_W1, decam_z

from IPython import embed

btsettl_path = get_project_root() / 'resource' / 'dwarf' / 'bt-settl'
setting = 'logg-5-metallicity-0'
files = [f for f in os.listdir(btsettl_path / setting) if f.startswith('lte')]
for f in files:
    dat = ascii.read(btsettl_path / setting / f)
    temp = int(f[3:6]) * 100
    wl, flux = dat['col1'], dat['col2']
    mask = (wl > 3000) & (wl < 60000)

    wl, flux = wl[mask] * u.AA, flux[mask] * u.erg / u.cm ** 2/ u.s / u.AA

    # write wl, flux to fits file
    dat = Table([wl, flux], names=['wave', 'flux'], units=['AA', 'erg / (cm2 s AA)'])
    dat.write(btsettl_path / setting / f'lte-{temp}K.fits', overwrite=True)