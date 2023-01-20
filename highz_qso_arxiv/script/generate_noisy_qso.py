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

sim_qso = pyfits.getdata('../resource/catalog/high_z_QSOs.fits', 1)
fz_simqso = sim_qso['DECam-DECaLS-z']
fW1_simqso = sim_qso['WISE-unWISE-W1']
fW2_simqso = sim_qso['WISE-unWISE-W2']
fY_simqso = sim_qso['VISTA-VISTA-Y']
fJ_simqso = sim_qso['VISTA-VISTA-J']
fH_simqso = sim_qso['VISTA-VISTA-H']
fK_simqso = sim_qso['VISTA-VISTA-K']

qso_flux = np.array([fJ_simqso, fz_simqso, fY_simqso, fH_simqso, fK_simqso, fW1_simqso, fW2_simqso])
qso_flux, qso_flux_err = add_noise(qso_flux)
np.save('../resource/catalog/qso_flux_noisy.npy', qso_flux)
np.save('../resource/catalog/qso_flux_err.npy', qso_flux_err)