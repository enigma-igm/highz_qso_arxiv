import os
import speclite
import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from scipy import interpolate
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval, ImageNormalize, SqrtStretch
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

sim_lris_path = '../resource/simulation/lris/'
targets = ['J1319+0101', 'J0901+2906', 'J1724+3718', 'J1326+0927']
trace_ids = [2, 2, 1, 2]
redshifts = [5.726, 6.093, 5.733, None]
m_Js = [21.20, 20.9, 21, 20.12]

i = 3
target = targets[i]
trace_id = trace_ids[i]
redshift = redshifts[i]
m_J = m_Js[i]

hdul = io.fits_open(os.path.join(sim_lris_path, f'{target}_tellmodel.fits'))
telluric = hdul[1].data

spec2dfile = os.path.join(sim_lris_path, f'spec2d_{target}.fits')
spec1dfile = spec2dfile.replace('spec2d', 'spec1d')

detname = DetectorContainer.get_name(det=1)
spec2DObj = spec2dobj.Spec2DObj.from_file(spec2dfile, detname, chk_version=False)
sobjs = specobjs.SpecObjs.from_fitsfile(spec1dfile, chk_version=False)
header = fits.getheader(spec1dfile)

gpm = spec2DObj.bpmmask == 0
image = (spec2DObj.sciimg - spec2DObj.skymodel) * np.sqrt(spec2DObj.ivarmodel) * gpm
mask = (image < 5) & (image > -5)
norm = ImageNormalize(image, interval=ZScaleInterval(), stretch=SqrtStretch())
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.imshow(image.T, origin='lower', norm=norm, cmap='gray')
ax.set_ylim(855, 1000)
ax.set_xticks([])
ax.set_yticks([])
plt.show()

