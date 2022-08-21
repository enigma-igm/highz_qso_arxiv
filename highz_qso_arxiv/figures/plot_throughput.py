import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from astropy.io import fits, ascii

from pypeit import io
from pypeit.sensfunc import IRSensFunc
from pypeit.core.flux_calib import Nlam_to_Flam
from highz_qso_arxiv import plot as hz_plot

from IPython import embed

fig, ax = plt.subplots(figsize=(10, 6))

# dat = ascii.read("../resource/sensfunc/keck_lris_throughput.dat")
# wave = dat["col1"]
# throughput = dat["col2"]

hdul = io.fits_open("../resource/sensfunc/GD153_lris_long_8.7_sens.fits")
sens = IRSensFunc.from_hdu(hdul)
wave = sens.wave
throughput = sens.throughput

ax.plot(wave/1e4, throughput, label="LRIS")

# dat = ascii.read("../resource/sensfunc/keck_mosfire_throughput.dat")
# wave = dat["col1"]
# throughput = dat["col2"]

hdul = io.fits_open("../resource/sensfunc/GD153_mosfire_sens.fits")
sens = IRSensFunc.from_hdu(hdul)
wave = sens.wave
throughput = sens.throughput
ax.plot(wave/1e4, throughput, label="MOSFIRE")

hdul = io.fits_open("../resource/sensfunc/GD153_nires_sens.fits")
sens = IRSensFunc.from_hdu(hdul)
wave = sens.wave
throughput = sens.throughput
for j in range(wave.shape[1]):
    _wave = wave[:, j]
    _throughput = throughput[:, j]
    mask = _throughput > 0
    label = f"NIRES order={7-j}"
    ax.plot(_wave[mask]/1e4, _throughput[mask], label=label)

ax.set_ylim(0, 0.5)
ax.set_xlabel(r'Wavelength [$\mu m$]', fontsize=20)
ax.set_ylabel('Throughput', fontsize=20)
ax.legend(fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=18)

wave_min, wave_max = ax.get_xlim()
wave_min, wave_max = wave_min*1e4, wave_max*1e4
z_min, z_max = wave_min / 1215.67 - 1, wave_max / 1215.67 - 1
ax2 = ax.twiny()
ax2.set_xlim(z_min, z_max)
ax2.set_xlabel("redshift", fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=18)

# qso_cat = pd.read_csv("../resource/machine_readable_quasar_list.1.15.csv")
# redshift = list(qso_cat[" z "])
# redshift.pop(redshift.index('>6.47'))
# redshift = [1215.67 * (1 + float(z)) for z in redshift]
# ax3 = ax.twinx()
# ax3.hist(redshift, color="black", bins=20, alpha=0.8, label=r"known $z>5.7$ quasars")
# ax3.set_ylim(1200, 0)
# ax3.vlines(1215.67 * (1 + 5.7), 0, 100, ls="dashed", color="black")
# ax3.hlines(100, 5500, 1215.67 * (1 + 5.7), ls="dashed", color="black")
# ax3.fill_between([5500, 1215.67 * (1 + 5.7)], [100, 100], [0, 0], color="black", alpha=0.5)
# ax3.set_xlim(5500, 25000)
# ax3.set_yticklabels([])
# ax3.legend(loc="lower right", fontsize=12)

# ax.set_xscale("log")
# ax2.set_xscale("log")
# ax3.set_xscale("log")
# plt.show()
plt.savefig("./throughput.pdf")