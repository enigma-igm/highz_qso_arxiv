import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
plt.style.use('science')
from pypeit.utils import calc_ivar

import warnings
warnings.filterwarnings("ignore")

def inverse(array):
    return (array > 0.0)/(np.abs(array) + (array == 0.0))
    
def ivarsmooth(flux, ivar, window):
    nflux = (flux.shape)[0]
    halfwindow = int(np.floor((np.round(window) - 1)/2))
    shiftarr = np.zeros((nflux, 2*halfwindow + 1))
    shiftivar = np.zeros((nflux, 2*halfwindow + 1))
    shiftindex = np.zeros((nflux, 2*halfwindow + 1))
    indexarr = np.arange(nflux)
    indnorm = np.outer(indexarr,(np.zeros(2 *halfwindow + 1) + 1))
    for i in np.arange(-halfwindow,halfwindow + 1,dtype=int):
        shiftarr[:,i+halfwindow] = np.roll(flux,i)
        shiftivar[:, i+halfwindow] = np.roll(ivar, i)
        shiftindex[:, i+halfwindow] = np.roll(indexarr, i)
    wh = (np.abs(shiftindex - indnorm) > (halfwindow+1))
    shiftivar[wh]=0.0
    outivar = np.sum(shiftivar,axis=1)
    nzero, = np.where(outivar > 0.0)
    zeroct=len(nzero)
    smoothflux = np.sum(shiftarr * shiftivar, axis=1)
    if(zeroct > 0):
        smoothflux[nzero] = smoothflux[nzero]/outivar[nzero]
    else:
        smoothflux = np.roll(flux, 2*halfwindow + 1)
    return (smoothflux, outivar)

def plot(name_list, fits_list, idx_list, ylim_list):
    num = len(fits_list)
    fig, axs = plt.subplots(num, 1, figsize=(12,3*num))
    if num == 1:
        axs = [axs]
    for idx, ax in enumerate(axs):
        hdul = fits.open(fits_list[idx])
        output = Table(hdul[idx_list[idx]].data)
        try:
            flux = np.array(output['OPT_FLAM'])
            flux_ivar = np.array(output['OPT_FLAM_IVAR'])
            flux_err = 1 / np.sqrt(flux_ivar)
            wave = np.array(output['OPT_WAVE'])
        except KeyError:
            flux = np.array(output['flux'])
            flux_ivar = np.array(output['ivar'])
            flux_err = 1 / np.sqrt(flux_ivar)
            wave = np.array(output['wave'])
            
        flux_sm, flux_ivar_sm = ivarsmooth(flux, flux_ivar, 7)

        ax.plot(wave[wave>5000], flux_sm[wave>5000], label=name_list[idx], color="black", lw=1.5)
        ax.plot(wave[wave>5000], calc_ivar(flux_ivar_sm)[wave>5000], 
                lw=1, color="red", alpha=0.6)
        ax.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
        ax.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
        # ax.set_title("J0847+0139, window=3", fontsize=20)
        ax.set_xlim(9900, 24000)
        ax.set_ylim(ylim_list[idx][0], ylim_list[idx][1])
        ax.legend(loc="upper right")
    fig.tight_layout()
    return fig, axs

"""
    LRIS
"""
fig, ax = plt.subplots(figsize=(20,6))

hdul = fits.open("J2212+2040_LRIS_coadd1d_tellcorr.fits")
output = Table(hdul[1].data)
wave, flux, flux_ivar = output["wave"], output["flux"], output["ivar"]
# mask = (wave<9400)
mask = wave>0
flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
ax.plot(wave[mask], flux_sm, color="black", alpha=0.8, lw=1, label="LRIS")
ax.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm)), lw=1, color="red", alpha=0.5)

hdul = fits.open("J2212+2040_LRIS_coadd1d_tellmodel.fits")
output = Table(hdul[1].data)
wave, flux = output["WAVE"].value[0], output["TELLURIC"].value[0]
ax2 = ax.twinx()
ax2.plot(wave[mask], flux[mask], color="navy", zorder=1, alpha=0.4)

ax2.set_ylim(-2,1.)
ax2.set_yticks([0,0.5,1])
ax.set_ylim(np.mean(flux_sm)-2*np.std(flux_sm),np.mean(flux_sm)+2*np.std(flux_sm))
# ax.set_xlim(7300, 9400)
ax.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
ax.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
plt.savefig("J2212+2040_LRIS.pdf")

"""
    MOSFIRE
"""
fig, ax = plt.subplots(figsize=(20,6))

hdul = fits.open("J2212_MOSFIRE_coadd_tellcorr.fits")
output = Table(hdul[1].data)
wave, flux, flux_ivar = output["wave"], output["flux"], output["ivar"]
mask = (wave>9700)
# mask = wave>0
flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
ax.plot(wave[mask], flux_sm, color="black", alpha=0.8, lw=1, label="MOSFIRE")
ax.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm)), lw=1, color="red", alpha=0.5)

hdul = fits.open("J2212_MOSFIRE_coadd_tellmodel.fits")
output = Table(hdul[1].data)
wave, flux = output["WAVE"].value[0], output["TELLURIC"].value[0]
ax2 = ax.twinx()
ax2.plot(wave[mask], flux[mask], color="navy", zorder=1, alpha=0.4)

ax2.set_ylim(-2,1.)
ax2.set_yticks([0,0.5,1])
ax.set_ylim(np.mean(flux_sm)-2*np.std(flux_sm),np.mean(flux_sm)+2*np.std(flux_sm))
# ax.set_xlim(7300, 9400)
ax.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
ax.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
plt.savefig("J2212+2040_MOSFIRE.pdf")

"""
    NIRES
"""
fig, ax = plt.subplots(figsize=(20,6))

hdul = fits.open("J2212+2040_NIRES_coadd1d_tellcorr.fits")
output = Table(hdul[1].data)
wave, flux, flux_ivar = output["wave"], output["flux"], output["ivar"]
# mask = (wave>10400)
mask = wave>0
flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
ax.plot(wave[mask], flux_sm, color="black", alpha=0.8, lw=1, label="NIRES")
ax.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm)), lw=1, color="red", alpha=0.5)

hdul = fits.open("J2212+2040_NIRES_coadd1d_tellmodel.fits")
output = Table(hdul[1].data)
wave, flux = output["WAVE"].value[0], output["TELLURIC"].value[0]
ax2 = ax.twinx()
ax2.plot(wave[mask], flux[mask], color="navy", zorder=1, alpha=0.4)

ax2.set_ylim(-2,1)
ax2.set_yticks([0,0.5,1])
ax.set_ylim(np.mean(flux_sm)-2*np.std(flux_sm),np.mean(flux_sm)+6*np.std(flux_sm))
ax.set_ylim(-0.05,0.8)
# ax.set_xlim(10400, 24500)
ax.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
ax.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
plt.savefig("J2212+2040_NIRES.pdf")

"""
    LRIS+MOSFIRE+NIRES
"""
fig, (ax2, ax3, ax1) = plt.subplots(3, 1, figsize=(20,20))

hdul = fits.open("J2212+2040_LRIS_coadd1d_tellcorr.fits")
output = Table(hdul[1].data)
wave = output["wave"]
flux = output["flux"]
flux_ivar = output["ivar"]
# mask = (wave<10400)
mask = wave>0
flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 1)
ax2.plot(wave[mask], flux_sm, color="black", lw=1, label="LRIS")
ax2.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm)), lw=1, color="red", alpha=0.5)
hdul = fits.open("J2212+2040_LRIS_coadd1d_tellmodel.fits")
output = Table(hdul[1].data)
wave = output["WAVE"].value[0]
flux = output["TELLURIC"].value[0]
ax2_twin = ax2.twinx()
ax2_twin.plot(wave[mask], flux[mask], color="grey", zorder=1, alpha=0.6)

hdul = fits.open("J2212_MOSFIRE_coadd_tellcorr.fits")
output = Table(hdul[1].data)
wave = output["wave"]
flux = output["flux"]
flux_ivar = output["ivar"]
# mask = (wave>9700)
mask = wave>0
flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
ax2.plot(wave[mask], flux_sm, color="slateblue", alpha=0.5, lw=1, label="MOSFIRE")

ax2.legend(loc="upper right")
ax2.set_ylim(-0.3,0.8)
ax2.set_xlim(8000, 14000)
ax2.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
ax2.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
ax2_twin.set_ylim(-2,1)
ax2_twin.set_yticks([0,0.5,1])

hdul = fits.open("J2212_MOSFIRE_coadd_tellcorr.fits")
output = Table(hdul[1].data)
wave = output["wave"]
flux = output["flux"]
flux_ivar = output["ivar"]
# mask = (wave>9700)
mask = wave>0
flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
ax3.plot(wave[mask], flux_sm, color="black", lw=1, label="MOSFIRE")
ax3.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm)), lw=1, color="red", alpha=0.5)
hdul = fits.open("J2212_MOSFIRE_coadd_tellmodel.fits")
output = Table(hdul[1].data)
wave = output["WAVE"].value[0]
flux = output["TELLURIC"].value[0]
ax3_twin = ax3.twinx()
ax3_twin.plot(wave[mask], flux[mask], color="grey", zorder=1, alpha=0.6)

ax3.legend(loc="upper right")
ax3.set_ylim(-0.3,0.8)
ax3.set_xlim(8000, 14000)
ax3.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
ax3.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
ax3_twin.set_ylim(-2,1)
ax3_twin.set_yticks([0,0.5,1])

hdul = fits.open("J2212+2040_NIRES_coadd1d_tellcorr.fits")
output = Table(hdul[1].data)
wave = output["wave"]
flux = output["flux"]
flux_ivar = output["ivar"]
# mask = ((wave>10400))
mask = wave>0
flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
ax1.plot(wave[mask], flux_sm, color="black", lw=1, label="NIRES")
ax1.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm)), lw=1, color="red", alpha=0.5)
hdul = fits.open("J2212+2040_NIRES_coadd1d_tellmodel.fits")
output = Table(hdul[1].data)
wave = output["WAVE"].value[0]
flux = output["TELLURIC"].value[0]
ax1_twin = ax1.twinx()
ax1_twin.plot(wave[mask], flux[mask], color="grey", zorder=1, alpha=0.6)

hdul = fits.open("J2212_MOSFIRE_coadd_tellcorr.fits")
output = Table(hdul[1].data)
wave = output["wave"]
flux = output["flux"]
flux_ivar = output["ivar"]
# mask = (wave>9700)
mask = wave>0
flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
ax1.plot(wave[mask], flux_sm, color="slateblue", alpha=0.5, lw=1, label="MOSFIRE")

ax1.legend(loc="upper left")
ax1.set_ylim(-0.3,0.8)
ax1.set_xlim(8000, 14000)
ax1.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
ax1.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
ax1_twin.set_ylim(-2,1)
ax1_twin.set_yticks([0,0.5,1])

plt.savefig("J2212_LRIS+MOSFIRE+NIRES.pdf")
