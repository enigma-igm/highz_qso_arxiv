
from cProfile import label
import os
import speclite
import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from scipy import interpolate
from astropy.stats import sigma_clipped_stats
from astropy.cosmology import FlatLambdaCDM, z_at_value

from highz_qso_arxiv.resource.filters import ukirt_J, hsc_z, decam_z, sdss_i
from highz_qso_arxiv.util.spec2dutil import gauss_comb
from highz_qso_arxiv.util import luminosity_to_flux, redshift_to_distance, get_project_root, ivarsmooth, rgb2gray
from highz_qso_arxiv.util.spec1dutil import add_gp_trough, add_damping_wing, add_telluric, extend_to_lower

from pypeit import io
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit.utils import inverse
from pypeit.core import extract
from pypeit.display import display
from pypeit.sensfunc import IRSensFunc
from pypeit.core.wavecal import wvutils
from pypeit.core.moment import moment1d
from pypeit.core.flux_calib import Flam_to_Nlam
from pypeit.spectrographs.util import load_spectrograph
from pypeit.images.detector_container import DetectorContainer
from pypeit.core import procimg

from IPython import embed

import matplotlib as mpl
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

root = get_project_root()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20,6), gridspec_kw={'height_ratios': [1, 4]})
trace = plt.imread(os.path.join(root, "resource", "MOSFIRE2204_J1332+0150.tiff"))
trace = rgb2gray(trace)
ax1.imshow(np.flip(trace.T), cmap='gray', aspect='auto')

ax1.set_xticklabels([])
ax1.set_yticklabels([])

spec0 = fits.getdata(os.path.join(root, "resource", "J1332+0150_coadd.fits"), 1)
waverange = (spec0['wave'][0], spec0['wave'][-1])
spec = fits.getdata(os.path.join(root, "resource", "spec1d_m220409-m220409-J1332+0150.fits"), 1)
mask = (spec['OPT_WAVE'] > waverange[0]) & (spec['OPT_WAVE'] < waverange[1])
smoothed_counts, smoothed_ivar = ivarsmooth(spec['OPT_COUNTS'][mask], spec['OPT_COUNTS_IVAR'][mask], 5)
ax2.plot(spec['OPT_WAVE'][mask], smoothed_counts, drawstyle='steps-mid', color=CB_color_cycle[0], lw=1.5)
ax2.plot(spec['OPT_WAVE'][mask], np.sqrt(inverse(smoothed_ivar)), drawstyle='steps-mid', color=CB_color_cycle[1], lw=1.5)


# ax2.plot(spec['OPT_WAVE'][mask], spec['OPT_FLAM'][mask], drawstyle='steps-mid', color=CB_color_cycle[2], lw=1.5)
# ax3 = ax2.twinx()
# ax3.plot(spec['OPT_WAVE'][mask], spec['OPT_COUNTS'][mask], drawstyle='steps-mid', color=CB_color_cycle[0], lw=1.5)
ax2.set_ylim(-5, 70)
# ax3.set_ylim(-10, 200)
ax2.set_xlabel('Wavelength [Angstrom]', fontsize=35)
ax2.set_ylabel('Flux', fontsize=35)
ax2.tick_params(labelsize=30)
ax2.legend(fontsize=30)
# ax2.text(0.18, 0.9,'J1319+0101 (QSO, z=5.726)', ha='center', va='center', transform=ax2.transAxes, fontsize=30)

plt.subplots_adjust(hspace=0.)
plt.show()