import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from scipy import interpolate

from pypeit.utils import inverse
from highz_qso_arxiv.util.spec1dutil import rescale
from highz_qso_arxiv.resource.filters import ukirt_J
from IPython import embed

hdul = fits.open('../arxiv/LRIS_2201/reduced/all/J1326+0927/J1326+0927_coadd.fits')
wl = hdul[1].data['wave']
flux = hdul[1].data['flux']
ivar = hdul[1].data['ivar']
err = np.sqrt(inverse(ivar))

temp_iras = ascii.read(f'../resource/dwarf/M9III_IRAS15060+0947.dat')
temp_lris = ascii.read(f'../resource/dwarf/keck_lris_M9_1.dat')
wl_iras, flux_iras = temp_iras["col1"] * 1e4, temp_iras["col2"]
wl_lris, flux_lris = temp_lris["col1"], temp_lris["col2"]

flux_iras = rescale(wl_iras, flux_iras, wl, flux, err)
flux_lris = rescale(wl_lris, flux_lris, wl, flux, err)

iras_mask = flux_iras > 0
flux_iras = flux_iras[iras_mask]
wl_iras = wl_iras[iras_mask]

# combine two spectrum
# the first part use lris spectrum
# when the wavelength is larger than lris, use iras spectrum
comb_mask = wl_iras > np.max(wl_lris)
wl_comb = np.append(wl_lris, wl_iras[comb_mask])
flux_comb = np.append(flux_lris, flux_iras[comb_mask])

# plot the template and the data
plt.plot(wl, flux, label="data", color='black')
plt.plot(wl, err, color='red')
# plt.plot(wl_iras, flux_iras, label="SpeX")
# plt.plot(wl_lris, flux_lris, label="LRIS")
plt.plot(wl_comb, flux_comb, label="combine")
plt.show()

print(ukirt_J.get_ab_magnitudes(flux_iras*1e-17*u.erg/u.s/u.cm**2/u.AA, wl_iras*u.AA))