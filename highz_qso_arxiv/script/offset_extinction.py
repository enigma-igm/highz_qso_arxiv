import imp
import pandas as pd
import numpy as np
from astropy.io import fits

from highz_qso_arxiv.util import mjd_to_unix
from highz_qso_arxiv.util.photutil import find_peak, get_mosfire_acq, get_mosfire_acq_proc, circle_mask, naive_bkg_subtract, sigma_clipping_bkg_subtract
from highz_qso_arxiv.crawler import get_skyprobe_extinction
from highz_qso_arxiv.resource.zero_points import mosfire_ZP
from highz_qso_arxiv.plot import plot_acq_and_hist, plot_extinction
from photutils.aperture import CircularAperture, RectangularAperture, aperture_photometry

from IPython import embed

path = "../arxiv/MOSFIRE_2204/raw"
df = pd.read_csv('offset_list.csv')
# df = df[df["flag"]==1]

m_aper = []
plate_scale = 0.1798
lengthx=round(2/plate_scale)
lengthy=round(3.5/plate_scale)
for i in range(len(df)):
    idx = df.index[i]
    acq_img = get_mosfire_acq_proc(path, objfile="m220409_{}.fits".format("%04d"%df["objframe"][idx]), 
                                   skyfile="m220409_{}.fits".format("%04d"%df["skyframe"][idx]))
    hdr = fits.getheader("{}/m220409_{}.fits".format(path, "%04d"%df["objframe"][idx]))
    peak = find_peak(acq_img)
    # plot_acq_and_hist(acq_img, peak, title=df["offset"][idx], display=True)

    mask_radius = 6. # pix
    star_medsub = naive_bkg_subtract(acq_img, mask_source=True, mask_radius=mask_radius)
    # star_sigmaclip = sigma_clipping_bkg_subtract(acq_img, mask_radius=10)

    # try use photutil here
    aperture = CircularAperture(peak, r=mask_radius)
    phot_table = aperture_photometry(star_medsub, aperture, method="exact")[0]

    # mask I use before using photutil
    mask = circle_mask(acq_img, peak[0], peak[1], mask_radius)

    f_star_photutil = phot_table["aperture_sum"] / hdr["TRUITIME"]
    f_star = np.sum(star_medsub[mask]) / hdr["TRUITIME"]
    # f_star = np.sum(star_sigmaclip[mask]) / hdr["TRUITIME"]

    m_star = -2.5 * np.log10(f_star) + mosfire_ZP["J"]
    m_star_photutil = -2.5 * np.log10(f_star_photutil) + mosfire_ZP["J"]
    m_aper.append(m_star_photutil)
    print(df["offset"][idx], df["JAB"][idx], m_star, m_star_photutil)
m_aper = np.array(m_aper)

extinction, scatter, time_mjd, time_unix, time_closest = [], [], [], [], []
for i in range(len(df)):
    idx = df.index[i]
    hdr = fits.getheader("{}/m220409_{}.fits".format(path, "%04d"%df["objframe"][idx]))
    obs_mjd = hdr["MJD-OBS"]
    obs_unix = mjd_to_unix(obs_mjd)

    time_mjd.append(obs_mjd)
    time_unix.append(obs_unix)
    
    extinction_data = get_skyprobe_extinction(obs_unix-900, obs_unix+900)
    closest_time = min(extinction_data["time"], key=lambda x:abs(x-obs_unix))
    dat = extinction_data[extinction_data["time"]==closest_time]
    extinction.append(dat["extinction"].value[0]-0.03)
    scatter.append(dat["scatter"].value[0])
    time_closest.append(closest_time)
extinction = np.array(extinction)
scatter = np.array(scatter)
time_mjd = np.array(time_mjd)
time_unix = np.array(time_unix)
time_closest = np.array(time_closest)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6,6))

flag_good = df["flag"]==1
flag_bad = df["flag"]==2
flag_feige = df["flag"]==0
ax.plot([0,20], [0,20])
ax.scatter(df["JAB"][flag_good], m_aper[flag_good], color="red", label=r"$\rm m_{aper}$")
# ax.scatter(df["JAB"][flag_good], m_aper-extinction[flag_good], label=r"$\rm m_{aper}-extinction$")

ax.scatter(df["JAB"][flag_bad], m_aper[flag_bad], color="grey", label="flag: near-edge")
# ax.scatter(df["JAB"][flag_bad], m_aper-extinction[flag_bad], label=r"$\rm m_{aper}-extinction$")

ax.scatter(df["JAB"][flag_feige], m_aper[flag_feige], color="blue", alpha=0.5, label="flag: feige")
# ax.scatter(df["JAB"][flag_feige], m_aper-extinction[flag_feige], label=r"$\rm m_{aper}-extinction$")


ax.set_xlim(14, 20)
ax.set_ylim(14, 20)
ax.set_xlabel(r"True $\rm m_{J}$", fontsize=15)
ax.set_ylabel(r"ACQ $\rm m_{J}$", fontsize=15)
ax.legend()
plt.show()

all_dat = get_skyprobe_extinction(time_unix[0]-900, time_unix[-1]+900)

fig, ax = plt.subplots(figsize=(10,6))
ax.scatter(time_unix[flag_good], m_aper[flag_good]-df["JAB"][flag_good], color="red", label=r"$m_{\rm aper}-m_{J}$")
ax.plot(time_unix[flag_good], m_aper[flag_good]-df["JAB"][flag_good], color="red")
ax.scatter(time_unix[flag_bad], m_aper[flag_bad]-df["JAB"][flag_bad], color="grey", label="flag: near-edge")
ax.scatter(time_unix[flag_feige], m_aper[flag_feige]-df["JAB"][flag_feige], color="blue", alpha=0.5, label="flag: feige")

ax.plot(all_dat["time"], all_dat["extinction"]-0.03, label="extinction")
ax.legend()
plt.show()