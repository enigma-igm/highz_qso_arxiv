from unicodedata import name
from highz_qso_arxiv.plot import plot_single, plot_series
from IPython import embed

plot_single("J2212 LRIS", "J2212+2040_LRIS_coadd1d_tellcorr.fits", 1, 
            telluric=True, telluric_fits_file="J2212+2040_LRIS_coadd1d_tellmodel.fits", 
            display=True, save_file="J2212+2040_LRIS.pdf")

plot_single("J2212 MOSFIRE", "J2212_MOSFIRE_coadd_tellcorr.fits", 1, 
            telluric=True, telluric_fits_file="J2212_MOSFIRE_coadd_tellmodel.fits", 
            display=True, save_file="J2212+2040_MOSFIRE.pdf")

plot_single("J2212 NIRES", "J2212+2040_NIRES_coadd1d_tellcorr.fits", 1, 
            telluric=True, telluric_fits_file="J2212+2040_NIRES_coadd1d_tellmodel.fits", 
            display=True, save_file="J2212+2040_NIRES.pdf")

# """
#     LRIS+MOSFIRE+NIRES
# """
# fig, (ax2, ax3, ax1) = plt.subplots(3, 1, figsize=(20,20))

# hdul = fits.open("J2212+2040_LRIS_coadd1d_tellcorr.fits")
# output = Table(hdul[1].data)
# wave = output["wave"]
# flux = output["flux"]
# flux_ivar = output["ivar"]
# # mask = (wave<10400)
# mask = wave>0
# flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 1)
# ax2.plot(wave[mask], flux_sm, color="black", lw=1, label="LRIS")
# ax2.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm)), lw=1, color="red", alpha=0.5)
# hdul = fits.open("J2212+2040_LRIS_coadd1d_tellmodel.fits")
# output = Table(hdul[1].data)
# wave = output["WAVE"].value[0]
# flux = output["TELLURIC"].value[0]
# ax2_twin = ax2.twinx()
# ax2_twin.plot(wave[mask], flux[mask], color="grey", zorder=1, alpha=0.6)

# hdul = fits.open("J2212_MOSFIRE_coadd_tellcorr.fits")
# output = Table(hdul[1].data)
# wave = output["wave"]
# flux = output["flux"]
# flux_ivar = output["ivar"]
# # mask = (wave>9700)
# mask = wave>0
# flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
# ax2.plot(wave[mask], flux_sm, color="slateblue", alpha=0.5, lw=1, label="MOSFIRE")

# ax2.legend(loc="upper right")
# ax2.set_ylim(-0.3,0.8)
# ax2.set_xlim(8000, 14000)
# ax2.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
# ax2.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
# ax2_twin.set_ylim(-2,1)
# ax2_twin.set_yticks([0,0.5,1])

# hdul = fits.open("J2212_MOSFIRE_coadd_tellcorr.fits")
# output = Table(hdul[1].data)
# wave = output["wave"]
# flux = output["flux"]
# flux_ivar = output["ivar"]
# # mask = (wave>9700)
# mask = wave>0
# flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
# ax3.plot(wave[mask], flux_sm, color="black", lw=1, label="MOSFIRE")
# ax3.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm)), lw=1, color="red", alpha=0.5)
# hdul = fits.open("J2212_MOSFIRE_coadd_tellmodel.fits")
# output = Table(hdul[1].data)
# wave = output["WAVE"].value[0]
# flux = output["TELLURIC"].value[0]
# ax3_twin = ax3.twinx()
# ax3_twin.plot(wave[mask], flux[mask], color="grey", zorder=1, alpha=0.6)

# ax3.legend(loc="upper right")
# ax3.set_ylim(-0.3,0.8)
# ax3.set_xlim(8000, 14000)
# ax3.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
# ax3.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
# ax3_twin.set_ylim(-2,1)
# ax3_twin.set_yticks([0,0.5,1])

# hdul = fits.open("J2212+2040_NIRES_coadd1d_tellcorr.fits")
# output = Table(hdul[1].data)
# wave = output["wave"]
# flux = output["flux"]
# flux_ivar = output["ivar"]
# # mask = ((wave>10400))
# mask = wave>0
# flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
# ax1.plot(wave[mask], flux_sm, color="black", lw=1, label="NIRES")
# ax1.plot(wave[mask], inverse(np.sqrt(flux_ivar_sm)), lw=1, color="red", alpha=0.5)
# hdul = fits.open("J2212+2040_NIRES_coadd1d_tellmodel.fits")
# output = Table(hdul[1].data)
# wave = output["WAVE"].value[0]
# flux = output["TELLURIC"].value[0]
# ax1_twin = ax1.twinx()
# ax1_twin.plot(wave[mask], flux[mask], color="grey", zorder=1, alpha=0.6)

# hdul = fits.open("J2212_MOSFIRE_coadd_tellcorr.fits")
# output = Table(hdul[1].data)
# wave = output["wave"]
# flux = output["flux"]
# flux_ivar = output["ivar"]
# # mask = (wave>9700)
# mask = wave>0
# flux_sm, flux_ivar_sm = ivarsmooth(flux[mask], flux_ivar[mask], 3)
# ax1.plot(wave[mask], flux_sm, color="slateblue", alpha=0.5, lw=1, label="MOSFIRE")

# ax1.legend(loc="upper left")
# ax1.set_ylim(-0.3,0.8)
# ax1.set_xlim(8000, 14000)
# ax1.set_xlabel(r"wavelength ($\AA$)", fontsize=15)
# ax1.set_ylabel(r"f$_{\lambda}$ ($10^{-17}$ ergs$^{-1}$cm$^{-2}\AA^{-1}$)", fontsize=15)
# ax1_twin.set_ylim(-2,1)
# ax1_twin.set_yticks([0,0.5,1])

# plt.savefig("J2212_LRIS+MOSFIRE+NIRES.pdf")
