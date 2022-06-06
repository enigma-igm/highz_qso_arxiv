import imp
import pandas as pd
import numpy as np
from astropy.io import fits

from highz_qso_arxiv.util import mjd_to_unix
from highz_qso_arxiv.util.photutil import find_peak, get_mosfire_acq, get_mosfire_acq_proc, circle_mask, naive_bkg_subtract, sigma_clipping_bkg_subtract
from highz_qso_arxiv.crawler import get_skyprobe_extinction
from highz_qso_arxiv.resource.zero_points import mosfire_ZP
from highz_qso_arxiv.plot import plot_acq_and_hist, plot_extinction

from IPython import embed

path = "../arxiv/MOSFIRE_2204/raw"
df = pd.read_csv('offset_list.csv')

m_aper = []
for i in range(len(df)):
    acq_img = get_mosfire_acq_proc(path, objfile="m220409_{}.fits".format("%04d"%df["objframe"][i]), 
                                   skyfile="m220409_{}.fits".format("%04d"%df["skyframe"][i]))
    hdr = fits.getheader("{}/m220409_{}.fits".format(path, "%04d"%df["objframe"][i]))
    peak = find_peak(acq_img)
    # plot_acq_and_hist(acq_img, peak, title=df["offset"][i], display=True)

    mask = circle_mask(acq_img, peak[0], peak[1], 6)
    star_medsub = naive_bkg_subtract(acq_img, mask_source=True)
    # star_sigmaclip = sigma_clipping_bkg_subtract(acq_img, mask_radius=10)

    f_star = np.sum(star_medsub[mask]) / hdr["TRUITIME"]
    # f_star = np.sum(star_sigmaclip[mask]) / hdr["TRUITIME"]

    m_star = -2.5 * np.log10(f_star) + mosfire_ZP["J"]
    m_aper.append(m_star)
    print(df["offset"][i], df["JAB"][i], m_star)
m_star = np.array(m_star)

extinction, scatter, time_mjd, time_unix = [], [], [], []
for i in range(len(df)):
    hdr = fits.getheader("{}/m220409_{}.fits".format(path, "%04d"%df["objframe"][i]))
    obs_mjd = hdr["MJD-OBS"]
    obs_unix = mjd_to_unix(obs_mjd)

    time_mjd.append(obs_mjd)
    time_unix.append(obs_unix)
    
    extinction_data = get_skyprobe_extinction(obs_unix-900, obs_unix+900)
    closest_time = min(extinction_data["time"], key=lambda x:abs(x-obs_unix))
    dat = extinction_data[extinction_data["time"]==closest_time]
    extinction.append(dat["extinction"].value[0]-0.03)
    scatter.append(dat["scatter"].value[0])
extinction = np.array(extinction)
scatter = np.array(scatter)
time_mjd = np.array(time_mjd)
time_unix = np.array(time_unix)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6,6))

ax.plot([0,20], [0,20])
ax.scatter(df["JAB"], m_aper-extinction)
ax.set_xlim(14, 20)
ax.set_ylim(14, 20)
ax.set_xlabel(r"$m_{J}$", fontsize=15)
ax.set_ylabel(r"$m_{\rm aper}-$extinction", fontsize=15)
ax.legend()
plt.show()

all_dat = get_skyprobe_extinction(time_unix[0]-900, time_unix[-1]+900)

fig, ax = plt.subplots(figsize=(10,6))
ax.scatter(time_unix, m_aper-df["JAB"], label=r"$m_{\rm aper}-m_{J}$")
ax.plot(time_unix, m_aper-df["JAB"])
ax.scatter(time_unix, extinction, label="extinction")
ax.plot(time_unix, extinction)
ax.plot(all_dat["time"], all_dat["extinction"], label="extinction")
ax.legend()
plt.show()