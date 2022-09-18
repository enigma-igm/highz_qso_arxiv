import corner
import numpy as np
import matplotlib as mpl
import mpl_scatter_density
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

from astropy.io import ascii
from astropy.table import Table
from matplotlib.patches import Ellipse

from highz_qso_arxiv import plot
from highz_qso_arxiv.plot import plot_contour
from highz_qso_arxiv.util import get_project_root, inverse
from highz_qso_arxiv.resource.filters import ukirt_J, wise_W1, decam_z

from IPython import embed

def error_cov(ax, xc, yc, cov, sigma=1,color='b', **kwargs):
    w, v = np.linalg.eigh(cov) # assumes symmetric matrix
    theta = np.degrees(np.arctan2(w[0]-cov[0,0],cov[0,1]))
    ellipse = Ellipse(xy=(xc,yc),width=2.*sigma*np.sqrt(w[0]),height=2.*sigma*np.sqrt(w[1]),angle=theta,color=color, zorder=2, **kwargs)
    ellipse.set_facecolor('none')
    ax.add_artist(ellipse)

uki_litqso = pyfits.getdata('../resource/catalog/qso_literature_uki.fits', 1)
redshift_litqso_uki = uki_litqso['redshift']
fz_litqso_uki, fz_err_litqso_uki = uki_litqso['flux_z_2'], uki_litqso['flux_err_z']
fW1_litqso_uki, fW1_err_litqso_uki = uki_litqso['flux_W1_2'], uki_litqso['flux_err_W1']
fJ_litqso_uki, fJ_err_litqso_uki = uki_litqso['flux_J'], uki_litqso['flux_J_err']

fig, ax = plt.subplots(1, 1, figsize=(8,8))
cov = np.ones((2, 2))
for i in range(len(fz_litqso_uki)):
    ax.scatter(fz_litqso_uki[i] / fJ_litqso_uki[i], fW1_litqso_uki[i] / fJ_litqso_uki[i], s=50, color='k')

    cov[0, 0] = (fz_litqso_uki[i] / fJ_litqso_uki[i])**2 * ((fz_err_litqso_uki[i] / fz_litqso_uki[i])**2 + (fJ_err_litqso_uki[i] / fJ_litqso_uki[i])**2)
    cov[1, 1] = (fW1_litqso_uki[i] / fJ_litqso_uki[i])**2 * ((fW1_err_litqso_uki[i] / fW1_litqso_uki[i])**2 + (fJ_err_litqso_uki[i] / fJ_litqso_uki[i])**2)
    cov[0, 1] = cov[1, 0] = fJ_err_litqso_uki[i]**2 / fJ_litqso_uki[i]**4
    error_cov(ax, fz_litqso_uki[i] / fJ_litqso_uki[i], fW1_litqso_uki[i] / fJ_litqso_uki[i], cov, sigma=1, color='k')
plt.show()

# cov = np.ones((2, 2))
# for i in range(len(relative_flux)):
#     for j in range(i + 1, len(relative_flux)):
#         ax = axes[j, i]
#         a = ax.scatter(relative_flux[i], relative_flux[j], marker='.', s=1, c=z, linewidths=0,
#                        cmap=discrete_cmap(N=5, base_cmap='autumn'), zorder=2)

#         for k in range(len(z)):
#             ax.scatter(relative_flux[i, k], relative_flux[j, k], marker=marker[k], s=60, linewidths=0, zorder=2,
#                        c=colors[k], label=label[k])
#             cov[0, 0] = np.square(relative_flux_err[i, k])
#             cov[1, 1] = np.square(relative_flux_err[j, k])
#             cov[1, 0] = cov[0, 1] = err_cov[k, i, j]
#             error_cov(ax=ax, xc=relative_flux[i, k], yc=relative_flux[j, k], cov=cov,
#                       color=colors[k])