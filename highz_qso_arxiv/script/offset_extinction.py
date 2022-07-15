import imp
import pandas as pd
import numpy as np
from astropy.io import fits

from highz_qso_arxiv.util import mjd_to_unix, ivarsmooth, inverse
from highz_qso_arxiv.util.photutil import find_peak, get_mosfire_acq, get_mosfire_acq_proc, circle_mask, naive_bkg_subtract, sigma_clipping_bkg_subtract
from highz_qso_arxiv.crawler import get_skyprobe_extinction
from highz_qso_arxiv.resource.zero_points import mosfire_ZP
from highz_qso_arxiv.plot import plot_acq_and_hist, plot_extinction
from photutils.aperture import CircularAperture, RectangularAperture, aperture_photometry

from IPython import embed

path = "../arxiv/MOSFIRE_2204/raw"
df = pd.read_csv('offset_star_list_2204.csv')
prefix = "m220409"
# prefix = "MF.20220111"
# df = df[df["flag"]==1]

m_aper = []
plate_scale = 0.1798
lengthx=round(2/plate_scale)
lengthy=round(3.5/plate_scale)

extinction, scatter, time_mjd, time_unix, time_closest = np.zeros(len(df)), np.zeros(len(df)), np.zeros(len(df)), np.zeros(len(df)), np.zeros(len(df))
m_aper, m_aper_err = np.zeros(len(df)), np.zeros(len(df))
airmass = np.zeros(len(df))
id = []
for i in range(len(df)):
    idx = df.index[i]

    """
        get aqcquisition image
    """
    # acq_img, acq_ivar = get_mosfire_acq_proc(path, objfile="{}.{}.fits.gz".format(prefix, "%04d"%df["objframe"][idx]), 
    #                                skyfile="{}.{}.fits.gz".format(prefix, "%04d"%df["skyframe"][idx]))
    acq_img, acq_ivar = get_mosfire_acq_proc(path, objfile="{}_{}.fits".format(prefix, "%04d"%df["objframe"][idx]), 
                                   skyfile="{}_{}.fits".format(prefix, "%04d"%df["skyframe"][idx]))
    acq_std = inverse(np.sqrt(acq_ivar))
    # hdr = fits.getheader("{}/{}.{}.fits.gz".format(path, prefix, "%04d"%df["objframe"][idx]))
    hdr = fits.getheader("{}/{}_{}.fits".format(path, prefix, "%04d"%df["objframe"][idx]))
    peak = find_peak(acq_img)
    # plot_acq_and_hist(acq_img, peak, title=df["offset"][idx], display=True)

    """
        do photometry
    """
    mask_radius = 6. # pix
    star_medsub = naive_bkg_subtract(acq_img, mask_source=True, mask_radius=mask_radius)
    # TODO: other background subtraction methods
    # star_sigmaclip = sigma_clipping_bkg_subtract(acq_img, mask_radius=10)

    aperture = CircularAperture(peak, r=mask_radius)
    phot_table = aperture_photometry(star_medsub, aperture, method="exact", error=acq_std)[0]

    f_star_photutil = phot_table["aperture_sum"] / hdr["TRUITIME"]
    m_star = -2.5 * np.log10(f_star_photutil) + mosfire_ZP["J"]
    m_star_err = np.abs(2.5*phot_table["aperture_sum_err"]/phot_table["aperture_sum"]/np.log(10))
    m_aper[i] = m_star
    m_aper_err[i] = m_star_err
    print(df["offset"][idx], df["jAperMag3"][idx], m_star, m_star-df["jAperMag3"][idx], m_star_err)
    id.append(df["offset"][idx])

    """
        get extinction data from SkyProbe
    """
    obs_mjd = hdr["MJD-OBS"]
    obs_unix = mjd_to_unix(obs_mjd)

    time_mjd[i] = obs_mjd
    time_unix[i] = obs_unix
    # search for extinction data +- 30min
    try:
        extinction_data = get_skyprobe_extinction(obs_unix-1800, obs_unix+1800)
    except ValueError:
        raise ValueError("{}: No extinction data found for {}-{}".
                         format(df["offset"][idx], obs_unix-1800, obs_unix+1800))

    closest_time = min(extinction_data["time"], key=lambda x:abs(x-obs_unix))
    dat = extinction_data[extinction_data["time"]==closest_time]
    extinction[i] = dat["extinction"].value[0]-0.03
    scatter[i] = dat["scatter"].value[0]
    time_closest[i] = closest_time

    airmass[i] = hdr["AIRMASS"]

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6,6))

flag_riccardo = df["flag"]==1 # riccardo and eduardo's targets
# flag_bad = df["flag"]==0 # at the edge or looks suspicious
flag_feige = df["flag"]==2
ax.plot([0,22], [0,22])
ax.scatter(df["jAperMag3"][flag_riccardo], m_aper[flag_riccardo], color="red", label=r"$\rm m_{aper}$")
ax.errorbar(df["jAperMag3"][flag_riccardo], m_aper[flag_riccardo], yerr=m_aper_err[flag_riccardo], color="red", fmt="o")
ax.scatter(df["jAperMag3"][flag_feige], m_aper[flag_feige], color="blue", label="flag: feige")
ax.errorbar(df["jAperMag3"][flag_feige], m_aper[flag_feige], yerr=m_aper_err[flag_feige], color="blue", fmt="o")
ax.set_xlim(np.min(df["jAperMag3"])-0.1, np.max(df["jAperMag3"])+0.1)
ax.set_ylim(np.min(m_aper)-0.1, np.max(m_aper)+0.1)
ax.set_xlabel(r"True $\rm m_{J}$", fontsize=15)
ax.set_ylabel(r"ACQ $\rm m_{J}$", fontsize=15)
ax.legend()
plt.show()

all_dat = get_skyprobe_extinction(time_unix[0]-2000, time_unix[-1]+2000)

fig, ax = plt.subplots(2, 1, figsize=(12,10))
ax[0].scatter(time_unix, m_aper-df["jAperMag3"], 
           color="red", label=r"$m_{\rm aper}-m_{J}$")
# put id on the plot
for i in range(len(id)):
    ax[0].text(time_unix[i], m_aper[i]-df["jAperMag3"][i], id[i], fontsize=10)
sorted_idx = np.argsort(time_unix)
ax[0].plot(time_unix[sorted_idx], m_aper[sorted_idx]-df["jAperMag3"][sorted_idx], color="grey")
# ax[0].scatter(time_unix[flag_feige], m_aper[flag_feige]-df["jAperMag3"][flag_feige], color="blue", label="flag: feige")

ax[1].plot(time_unix[sorted_idx], airmass[sorted_idx], color="grey")

ax[0].set_xlabel("UNIX TIME", fontsize=15)
ax[0].set_ylabel(r"$m_{\rm aper}-m_{J}$", fontsize=15)
ax[1].set_ylabel(r"$\rm Airmass$", fontsize=15)
ax[1].set_xlabel("UNIX TIME", fontsize=15)

ax[0].plot(all_dat["time"], all_dat["extinction"]-0.03, label="extinction")
ax[0].legend()
fig.tight_layout()
plt.show()