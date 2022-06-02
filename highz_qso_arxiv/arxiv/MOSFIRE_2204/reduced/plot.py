import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
plt.style.use('science')

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

hdul = fits.open("J2/Science_coadd/spec1d_m220409-m220409-J0841+3814.fits")
output = Table(hdul[7].data)
flux = np.array(output['OPT_FLAM'])
flux_ivar = np.array(output['OPT_FLAM_IVAR'])
flux_err = 1 / np.sqrt(flux_ivar)
wave = np.array(output['OPT_WAVE'])
flux_sm, flux_ivar_sm = ivarsmooth(flux, flux_ivar, 3)

fig, ax = plt.subplots(figsize=(10,5))
ax.plot(wave[wave>11000], flux_sm[wave>11000], label=r"smoothed $\rm OPT\_FLAM$")
ax.plot(wave[wave>11000], (1/np.sqrt(flux_ivar_sm)-np.min(1/np.sqrt(flux_ivar_sm))+0.05)[wave>11000], lw=0.5, color="grey")
ax.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
ax.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
ax.set_title("J0841+3814, window=3", fontsize=20)
ax.set_xlim(11200, 12500)
ax.set_ylim(0, 0.9)
ax.legend()
plt.savefig("J0841.pdf")