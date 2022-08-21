import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.table import Table
from astropy.io import ascii, fits

from pypeit.utils import inverse
from highz_qso_arxiv.util.spec1dutil import rescale
from highz_qso_arxiv.resource.filters import ukirt_J
from IPython import embed

for tp in ['M4.5', 'M5', 'M6', 'M7', 'M8', 'M9', 'M9.5', 'L0.5', 'L1', 'L2', 'L3', 'L5']:
    temp_irtf = ascii.read(f'irtf_{tp}_1.dat')
    temp_lris = ascii.read(f'keck_lris_{tp}_1.dat')
    wl_irtf, flux_irtf, err_irtf = temp_irtf["col1"] * 1e4, temp_irtf["col2"], temp_irtf["col3"]
    wl_lris, flux_lris = temp_lris["col1"], temp_lris["col2"]

    flux_lris, _ = rescale(wl_lris, flux_lris, wl_irtf, flux_irtf, err_irtf)

    irtf_mask = flux_irtf > 0
    flux_irtf = flux_irtf[irtf_mask]
    wl_irtf = wl_irtf[irtf_mask]

    # combine two spectrum

    # the first part use lris spectrum
    # when the wavelength is larger than lris, use irtf spectrum
    # use this in simulation
    comb_mask_lris = wl_irtf > np.max(wl_lris)
    wl_comb_lris = np.append(wl_lris, wl_irtf[comb_mask_lris])
    flux_comb_lris = np.append(flux_lris, flux_irtf[comb_mask_lris])

    # another combined way
    # when the wavelength is smaller than irtf, use lris spectrum
    # use this to measure photometry (zJW1)
    comb_mask_irtf = wl_lris < np.min(wl_irtf)
    wl_comb_irtf = np.append(wl_lris[comb_mask_irtf], wl_irtf)
    flux_comb_irtf = np.append(flux_lris[comb_mask_irtf], flux_irtf)

    # plot the template and the data
    # plt.plot(wl_comb_lris, flux_comb_lris, label="combine lris")
    plt.plot(wl_comb_irtf, flux_comb_irtf, label="combine irtf")
    plt.plot(wl_irtf, flux_irtf, 'r-', alpha=0.5, label='irtf')
    plt.show()

    dat = Table([wl_comb_lris, flux_comb_lris], names=['wave', 'flux'], units=['angstrom', 'W/m^2/um'])
    dat.write(f'combine_lris_{tp}.dat', format='ascii.ecsv', overwrite=True)

    dat = Table([wl_comb_lris, flux_comb_lris], names=['wave', 'flux'], units=['angstrom', 'W/m^2/um'])
    dat.write(f'combine_irtf_{tp}.dat', format='ascii.ecsv', overwrite=True)

