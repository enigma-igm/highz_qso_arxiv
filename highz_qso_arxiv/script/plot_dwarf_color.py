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

# dwarf_path = get_project_root() / 'resource' / 'dwarf'
# Mdwarf_type = ['M4.5', 'M5', 'M6', 'M7', 'M8', 'M9', 'M9.5', 'L0.5', 'L1', 'L2', 'L3', 'L5']
# czJ, cJW1 = [], []
# for tp in Mdwarf_type:
#     spec_irtf = ascii.read(dwarf_path / f'combine_irtf_{tp}.dat')
#     wl_irtf, flux_irtf = spec_irtf['wave'], spec_irtf['flux']
#     flux, wl = decam_z.pad_spectrum(flux_irtf, wl_irtf)
#     mask = flux < 0
#     flux, wl = flux[~mask], wl[~mask]
#     zJ = decam_z.get_ab_magnitudes(flux, wl)[0][0] - ukirt_J.get_ab_magnitudes(flux, wl)[0][0]
#     JW1 = ukirt_J.get_ab_magnitudes(flux, wl)[0][0] - wise_W1.get_ab_magnitudes(flux, wl)[0][0]
#     czJ.append(zJ)
#     cJW1.append(JW1)

# plt.scatter(czJ, cJW1, label='IRTF template')
# plt.plot(czJ, cJW1, ls='dashed', alpha=0.5)
# plt.xlabel('z-J', fontsize=25)
# plt.ylabel('J-W1', fontsize=25)
# # label Mdwarf_type along each point
# for i, tp in enumerate(Mdwarf_type):
#     plt.annotate(str(tp), (czJ[i], cJW1[i]), fontsize=15)

# plt.show()

czJ, cJW1 = [], []
temp = []
btsettl_path = get_project_root() / 'resource' / 'dwarf' / 'bt-settl'
setting = 'logg-3-metallicity-0'
files = [f for f in os.listdir(btsettl_path / setting) if f.startswith('lte-')]
for f in files:
    dat = pyfits.getdata(btsettl_path / setting / f, 1)
    wl, flux = dat['wave'], dat['flux']
    wl, ind = np.unique(wl, return_index=True)
    flux = flux[ind]
    
    # func = interp1d(wl, flux, kind='cubic')
    # wl = np.arange(3200, 55000, 1)
    # flux = func(wl)

    wl, flux = wl * u.AA, flux * u.erg / u.cm ** 2/ u.s / u.AA

    try:
        zJ = decam_z.get_ab_magnitudes(flux, wl)[0][0] - ukirt_J.get_ab_magnitudes(flux, wl)[0][0]
        JW1 = ukirt_J.get_ab_magnitudes(flux, wl)[0][0] - wise_W1.get_ab_magnitudes(flux, wl)[0][0]
    except ValueError:
        embed()
    czJ.append(zJ)
    cJW1.append(JW1)
    # extract the number between 'lte-' and 'K.fits'
    temp.append(int(f[4:-8]))

index = np.argsort(temp)
czJ = np.array(czJ)[index]
cJW1 = np.array(cJW1)[index]

plt.scatter(czJ, cJW1, label='BT-Settl')
plt.plot(czJ, cJW1, ls='dashed', alpha=0.5)
# add temp as text on the plot
for i, t in enumerate(temp):
    plt.annotate(str(t), (czJ[i], cJW1[i]), fontsize=15)

plt.xlabel('z-J', fontsize=25)
plt.ylabel('J-W1', fontsize=25)

plt.show()
