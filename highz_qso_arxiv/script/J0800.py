from highz_qso_arxiv.plot import plot_single
import matplotlib.pyplot as plt

# fits_file = '..//arxiv/NIRES_1903/reduced/all/coadd2d/J0800+3034_coadd_tellcorr.fits'
fits_file = '..//arxiv/MOSFIRE_2010/reduced/all/coadd2d/J0220+3458_coadd_tellcorr.fits'

fig, ax = plot_single('J0800', idx=1, fits_file=fits_file, 
                      smooth_window=11, display=False, template=False, telluric=True)
lines = [1216, 1549, 1909, 2799, 3727, 4862, 5008, 6564]
redshift = 7
for line in lines:
    ax.axvline(line * (1 + redshift), color='k', linestyle='--', alpha=1)
ax.set_ylim(-0.2, 1)
plt.show()