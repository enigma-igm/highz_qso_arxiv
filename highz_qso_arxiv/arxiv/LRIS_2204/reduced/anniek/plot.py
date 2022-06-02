import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
plt.style.use('science')

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
    fig, axs = plt.subplots(num, 1, figsize=(10,3*num))
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
            
        flux_sm, flux_ivar_sm = ivarsmooth(flux, flux_ivar, 3)

        ax.plot(wave[wave>5000], flux_sm[wave>5000], label=name_list[idx], color="black", lw=1.5)
        ax.plot(wave[wave>5000], inverse(flux_ivar_sm)[wave>5000], 
                lw=1, color="red", alpha=0.6)
        ax.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
        ax.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
        # ax.set_title("J0847+0139, window=3", fontsize=20)
        ax.set_xlim(7300,9600)
        ax.set_ylim(ylim_list[idx][0], ylim_list[idx][1])
        ax.legend(loc="upper left")
    fig.tight_layout()
    return fig, axs

name_list = ["J1401+4542", "J1523+2935", "J1535+6146"]
fits_list = ["reduced/spec1d_r220423_00100-J1401+4542_OFF_LRISr_20220423T131547.837.fits",
             "reduced/spec1d_r220423_00084-J1523+2935_OFF_LRISr_20220423T113251.619.fits",
             "reduced/J1535+6146_coadd.fits"]
idx_list = [2, 1, 1]
ylim_list = [(-0.1,0.8), (-0.1,2.5), (-0.05,0.38)]
fig, axs = plot(name_list, fits_list, idx_list, ylim_list)
fig.suptitle("LRIS 2204 (Anniek)", fontsize=25)
fig.subplots_adjust(top=0.9)
# plt.show()
plt.savefig("LRIS_2204_anniek.pdf")