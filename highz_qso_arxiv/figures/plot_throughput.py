import os
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
from scipy import interpolate

from pypeit import io
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit.utils import inverse
from pypeit.core import extract
from pypeit.display import display
from pypeit.sensfunc import IRSensFunc
from pypeit.core.wavecal import wvutils
from pypeit.core.moment import moment1d
from pypeit.core.flux_calib import Flam_to_Nlam
from pypeit.spectrographs.util import load_spectrograph
from pypeit.images.detector_container import DetectorContainer

from IPython import embed

import matplotlib as mpl

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
CB_color_cycle = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3']

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 

import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--plot', default=False, action='store_true')
parser.add_argument('--show', default=False, action='store_true')

args = parser.parse_args()

def simple_slit_loss(slitwidth, fwhm, platescale, binning):
    from scipy.special import erf
    return erf((slitwidth/2)/(np.sqrt(2)*fwhm*platescale*binning/2.355))

sim_path = '../resource/simulation'
spec1dfile = os.path.join(sim_path, 'nires', f'spec1d_GD153.fits')
sobjs_nires = specobjs.SpecObjs.from_fitsfile(spec1dfile, chk_version=False)
spec1dfile = os.path.join(sim_path, 'mosfire', f'spec1d_GD153.fits')
sobjs_mosfire = specobjs.SpecObjs.from_fitsfile(spec1dfile, chk_version=False)
spec1dfile = os.path.join(sim_path, 'lris', f'spec1d_GD153.fits')
sobjs_lris = specobjs.SpecObjs.from_fitsfile(spec1dfile, chk_version=False)


fig, ax = plt.subplots(figsize=(10, 8))

# hdul = io.fits_open("../resource/sensfunc/GD153_lris_long_8.7_sens.fits")
# sens = IRSensFunc.from_hdu(hdul)
# wave = sens.wave
# throughput = sens.throughput
# mask = throughput > 0
# ax.plot(wave[mask]/1e4, throughput[mask], color='red', lw=2, alpha=0.8, label=r"$\textbf{LRIS Red}$")

hdul = io.fits_open("../resource/sensfunc/GD153_lris_sens.fits")
sens = IRSensFunc.from_hdu(hdul)
wave = sens.wave
throughput = sens.throughput
slit_trans = simple_slit_loss(1, sobjs_lris[0].FWHMFIT, 0.123, 2)
seeing = 0.9
scaling = seeing / (np.mean(sobjs_lris[0].FWHMFIT) * 0.123 * 2)
pseudo_fwhm = sobjs_lris[0].FWHMFIT * scaling
slit_trans_psf = simple_slit_loss(1, pseudo_fwhm, 0.123, 2)

func_slit_trans = interpolate.interp1d(sobjs_lris[0].BOX_WAVE, slit_trans, fill_value="extrapolate")
func_slit_trans_psf = interpolate.interp1d(sobjs_lris[0].BOX_WAVE, slit_trans_psf, fill_value="extrapolate")

slit_loss = func_slit_trans(wave)
slit_loss_psf = func_slit_trans_psf(wave)

mask = throughput > 0
# ax.plot(wave[mask]/1e4, throughput[mask], color=CB_color_cycle[0], lw=2, alpha=0.8, label=r"$\textbf{LRIS Red}$")
# ax.plot(wave[mask]/1e4, throughput[mask] / slit_loss[mask] * slit_loss_psf[mask], color=CB_color_cycle[0], lw=2., alpha=0.8, ls="dashed")
ax.plot(wave[mask]/1e4, throughput[mask] / slit_loss[mask] * slit_loss_psf[mask], color=CB_color_cycle[0], lw=2., label=r"$\textbf{LRIS Red (seeing=0.9'')}$")
ax.plot(wave[mask]/1e4, throughput[mask] / slit_loss[mask], color=CB_color_cycle[0], lw=2., ls='dotted')

# dat = ascii.read("../resource/sensfunc/keck_mosfire_throughput.dat")
# wave = dat["col1"]
# throughput = dat["col2"]

hdul = io.fits_open("../resource/sensfunc/GD153_mosfire_sens.fits")
sens = IRSensFunc.from_hdu(hdul)
wave = sens.wave
throughput = sens.throughput
slit_trans = simple_slit_loss(1, sobjs_mosfire[0].FWHMFIT, 0.1798, 1)
seeing = 0.9
scaling = seeing / (np.mean(sobjs_mosfire[0].FWHMFIT) * 0.1798 * 1)
pseudo_fwhm = sobjs_mosfire[0].FWHMFIT * scaling
slit_trans_psf = simple_slit_loss(1, pseudo_fwhm, 0.1798, 1)

func_slit_trans = interpolate.interp1d(sobjs_mosfire[0].BOX_WAVE, slit_trans, fill_value="extrapolate")
func_slit_trans_psf = interpolate.interp1d(sobjs_mosfire[0].BOX_WAVE, slit_trans_psf, fill_value="extrapolate")

slit_loss = func_slit_trans(wave)
slit_loss_psf = func_slit_trans_psf(wave)

# ax.plot(wave/1e4, throughput, color=CB_color_cycle[1], lw=2, alpha=0.8, label=r"$\textbf{MOSFIRE-Y}$")
ax.plot(wave/1e4, throughput / slit_loss * slit_loss_psf, color=CB_color_cycle[1], lw=2., label=r"$\textbf{MOSFIRE-Y (seeing=0.9'')}$")
# ax.plot(wave/1e4, throughput / slit_loss * slit_loss_psf, color=CB_color_cycle[1], lw=2., alpha=0.8, ls="dashed")
ax.plot(wave/1e4, throughput / slit_loss, color=CB_color_cycle[1], lw=2., ls="dotted")

hdul = io.fits_open("../resource/sensfunc/GD153_nires_sens.fits")
sens = IRSensFunc.from_hdu(hdul)
wave = sens.wave
throughput = sens.throughput
for j in range(wave.shape[1]):
    _wave = wave[:, j]
    _throughput = throughput[:, j]
    slit_trans = simple_slit_loss(0.55, sobjs_nires[j*3].FWHMFIT, 0.15, 1)
    func_slit_trans = interpolate.interp1d(sobjs_nires[j*3].BOX_WAVE, slit_trans, fill_value="extrapolate")
    _slit_loss = func_slit_trans(_wave)

    mask = _throughput > 0
    if j < wave.shape[1] - 1:
        ax.plot(_wave[mask]/1e4, _throughput[mask], lw=2, color=CB_color_cycle[2], alpha=0.8)
        ax.plot(_wave[mask]/1e4, _throughput[mask] / _slit_loss[mask], lw=2., color=CB_color_cycle[2], ls="dotted")
        ax.text(_wave[mask][-1]/1e4, _throughput[mask][-1] / _slit_loss[mask][-1], f"{7-j}", color=CB_color_cycle[2], fontsize=15)
    else:
        lb = r"$\textbf{NIRES (order 3-7)}$" + "\n" + r"$\textbf{(order 7 seeing=0.9'')}$"
        ax.plot(_wave[mask]/1e4, _throughput[mask], lw=2, color=CB_color_cycle[2], alpha=0.8, label=lb)
        ax.plot(_wave[mask]/1e4, _throughput[mask] / _slit_loss[mask], lw=2., color=CB_color_cycle[2], ls="dotted")
        ax.text(_wave[mask][-1]/1e4, _throughput[mask][-1] / _slit_loss[mask][-1], f"{7-j}", color=CB_color_cycle[2], fontsize=15)

ax.legend(fontsize=16, loc="upper center")

ax_dummy = ax.twinx()
# ax_dummy.plot([], [], color="grey", lw=2., ls="-", label="Seeing = 0.9 arcsec \n (LRIS, MOSFIRE, NIRES order 7)")
ax_dummy.plot([], [], color="grey", lw=2., ls="dotted", label="Slit loss corrected \n assuming perfect centering")
ax_dummy.legend(fontsize=20, loc="lower right")
ax_dummy.set_yticks([])

ax.set_ylim(-0.01, 0.3)
ax.set_xlabel(r'Wavelength ($\mu m$)', fontsize=20)
ax.set_ylabel('Throughput', fontsize=20)
# ax.legend(fontsize=18, loc="lower right")
ax.tick_params(axis='both', which='major', labelsize=20, width=1, size=6)
ax.tick_params(axis='both', which='minor', labelsize=20, width=1, size=3)


wave_min, wave_max = ax.get_xlim()
wave_min, wave_max = wave_min*1e4, wave_max*1e4
z_min, z_max = wave_min / 1215.67 - 1, wave_max / 1215.67 - 1
ax2 = ax.twiny()
ax2.set_xlim(z_min, z_max)
ax2.set_xlabel(r"Ly$\alpha$ redshift", fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20, width=1, size=6)
ax2.tick_params(axis='both', which='minor', labelsize=20, width=1, size=3)

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

if args.show:
    plt.show()
else:
    plt.savefig("./throughput.pdf")