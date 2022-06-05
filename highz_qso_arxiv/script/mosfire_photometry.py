import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from astropy.io import fits

from highz_qso_arxiv.util import mjd_to_unix
from highz_qso_arxiv.util.photutil import find_peak, r_theta, get_mosfire_acq, circle_mask, naive_bkg_subtract, sigma_clipping_bkg_subtract
from highz_qso_arxiv.crawler import get_skyprobe_extinction
from highz_qso_arxiv.resource.zero_points import mosfire_ZP
from highz_qso_arxiv.plot import plot_acq, plot_hist, plot_acq_and_hist, plot_extinction

path = "../resource"
# 204 - star sky / 205 - star obj
# 206 - qso sky / 207 - qso obj
star_acq = get_mosfire_acq(path, "m220111_0204.fits", "m220111_0205.fits")
qso_acq = get_mosfire_acq(path, "m220111_0206.fits", "m220111_0207.fits")
qso_hdr = fits.getheader(f"{path}/m220111_0207.fits")
star_hdr = fits.getheader(f"{path}/m220111_0205.fits")

star_peak = find_peak(star_acq)
qso_peak = find_peak(qso_acq)

plot_acq_and_hist(star_acq, star_peak, display=True)
plot_acq_and_hist(qso_acq, qso_peak, display=True)

star_medsub = naive_bkg_subtract(star_acq, method="biweight")
star_sigmaclip = sigma_clipping_bkg_subtract(star_acq, mask_radius=10)

mask = circle_mask(star_medsub, star_peak[0], star_peak[1], 6)
f_star = np.sum(star_medsub[mask]) / star_hdr["ELAPTIME"]
m_star = -2.5 * np.log10(f_star) + mosfire_ZP["J"]
print(f"Star Magnitute with median subtraction: {m_star}")
f_star = np.sum(star_sigmaclip[mask]) / star_hdr["ELAPTIME"]
m_star = -2.5 * np.log10(f_star) + mosfire_ZP["J"]
print(f"Star Magnitute with sigma clipping: {m_star}")

qso_medsub = naive_bkg_subtract(qso_acq, method="biweight")
qso_sigmaclip = sigma_clipping_bkg_subtract(qso_acq, mask_radius=10)

mask = circle_mask(qso_medsub, qso_peak[0], qso_peak[1], 6)
f_qso = np.sum(qso_medsub[mask]) / qso_hdr["ELAPTIME"]
m_qso = -2.5 * np.log10(f_qso) + mosfire_ZP["J"]
print(f"qso Magnitute with median subtraction: {m_qso}")
f_qso = np.sum(qso_sigmaclip[mask]) / qso_hdr["ELAPTIME"]
m_qso = -2.5 * np.log10(f_qso) + mosfire_ZP["J"]
print(f"qso Magnitute with sigma clipping: {m_qso}")

qso_mjd = qso_hdr["MJD-obs"]
star_mjd = star_hdr["MJD-obs"]

qso_unix = mjd_to_unix(qso_mjd)
star_unix = mjd_to_unix(star_mjd)

extinction_data = get_skyprobe_extinction(star_unix-300, qso_unix+300)
print(f"star observing time: {star_unix}")
print(f"qso observing time: {qso_unix}")
print(extinction_data)

plot_extinction(extinction_data, star_unix, qso_unix, offset=0.03, display=True)