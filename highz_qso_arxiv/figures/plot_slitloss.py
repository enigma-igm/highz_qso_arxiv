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
from highz_qso_arxiv.plot import args
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

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=CB_color_cycle) 

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


fig, ax = plt.subplots(figsize=(10, 6))

# dat = ascii.read("../resource/sensfunc/keck_lris_throughput.dat")
# wave = dat["col1"]
# throughput = dat["col2"]

# hdul = io.fits_open("../resource/sensfunc/GD153_lris_long_8.7_sens.fits")
hdul = io.fits_open("../resource/sensfunc/GD153_lris_sens.fits")
sens = IRSensFunc.from_hdu(hdul)
wave = sens.wave
throughput = sens.throughput
slit_trans = simple_slit_loss(1, sobjs_lris[0].FWHMFIT, 0.123, 2)
seeing = 0.9
pseudo_fwhm = np.ones_like(sobjs_lris[0].FWHMFIT) * seeing / 0.123 / 2
slit_trans_psf = simple_slit_loss(1, pseudo_fwhm, 0.123, 2)

func_slit_trans = interpolate.interp1d(sobjs_lris[0].BOX_WAVE, slit_trans, fill_value="extrapolate")
func_slit_trans_psf = interpolate.interp1d(sobjs_lris[0].BOX_WAVE, slit_trans_psf, fill_value="extrapolate")

slit_loss = func_slit_trans(wave)
slit_loss_psf = func_slit_trans_psf(wave)

mask = throughput > 0
# ax.plot(wave[mask]/1e4, throughput[mask], color=CB_color_cycle[0], lw=2, alpha=0.5, label=r"$\textbf{LRIS Red}$")
# ax.plot(wave[mask]/1e4, throughput[mask] / slit_loss[mask] * slit_loss_psf[mask], color=CB_color_cycle[0], lw=1, alpha=0.5, ls="dashed")
# ax.plot(wave[mask]/1e4, throughput[mask] / slit_loss[mask], color=CB_color_cycle[0], lw=1, ls='-.')

# dat = ascii.read("../resource/sensfunc/keck_mosfire_throughput.dat")
# wave = dat["col1"]
# throughput = dat["col2"]

hdul = io.fits_open("../resource/sensfunc/GD153_mosfire_sens.fits")
sens = IRSensFunc.from_hdu(hdul)
wave = sens.wave
throughput = sens.throughput
slit_trans = simple_slit_loss(1, sobjs_mosfire[0].FWHMFIT, 0.1798, 1)
seeing = 0.9
pseudo_fwhm = np.ones_like(sobjs_mosfire[0].FWHMFIT) * seeing / 0.1798 / 1
slit_trans_psf = simple_slit_loss(1, pseudo_fwhm, 0.1798, 1)

func_slit_trans = interpolate.interp1d(sobjs_mosfire[0].BOX_WAVE, slit_trans, fill_value="extrapolate")
func_slit_trans_psf = interpolate.interp1d(sobjs_mosfire[0].BOX_WAVE, slit_trans_psf, fill_value="extrapolate")

slit_loss = func_slit_trans(wave)
slit_loss_psf = func_slit_trans_psf(wave)

# ax.plot(wave/1e4, throughput, color=CB_color_cycle[1], lw=2, alpha=0.5, label=r"$\textbf{MOSFIRE-Y}$")
# ax.plot(wave/1e4, throughput / slit_loss * slit_loss_psf, color=CB_color_cycle[1], lw=1, alpha=0.5, ls="dashed")
# ax.plot(wave/1e4, throughput / slit_loss, color=CB_color_cycle[1], lw=1, ls="-.")

hdul = io.fits_open("../resource/sensfunc/GD153_nires_sens.fits")
sens = IRSensFunc.from_hdu(hdul)
wave = sens.wave
throughput = sens.throughput

fhwm_onemicron = sobjs_nires[0].FWHMFIT[(np.abs(sobjs_nires[0].BOX_WAVE - 1e4) < 1)][0]
fhwm_onemicron = fhwm_onemicron * 0.15
seeing = 0.6
scaling = seeing / fhwm_onemicron
for j in range(wave.shape[1]):
    _wave = wave[:, j]
    _throughput = throughput[:, j]
    slit_trans = simple_slit_loss(0.55, sobjs_nires[j*3].FWHMFIT, 0.15, 1)
    pseudo_fwhm = sobjs_nires[j*3].FWHMFIT * scaling
    # pseudo_fwhm = np.ones_like(sobjs_nires[j*3].FWHMFIT) * seeing / 0.15 / 1
    slit_trans_psf = simple_slit_loss(0.55, pseudo_fwhm, 0.15, 1)

    func_slit_trans = interpolate.interp1d(sobjs_nires[j*3].BOX_WAVE, slit_trans, fill_value="extrapolate")
    func_slit_trans_psf = interpolate.interp1d(sobjs_nires[j*3].BOX_WAVE, slit_trans_psf, fill_value="extrapolate")
    _slit_loss = func_slit_trans(_wave)
    _slit_loss_psf = func_slit_trans_psf(_wave)

    mask = _throughput > 0
    if j < wave.shape[1] - 1:
        ax.plot(_wave[mask]/1e4, _throughput[mask], lw=2, color=CB_color_cycle[2], alpha=0.5)
        ax.plot(_wave[mask]/1e4, _throughput[mask] / _slit_loss[mask], lw=1, color='black', ls="-.")
        ax.plot(_wave[mask]/1e4, _throughput[mask] / _slit_loss[mask] * _slit_loss_psf[mask], lw=1, color='black', ls="dotted")
        ax.text(_wave[mask][-1]/1e4, _throughput[mask][-1], f"{7-j}", color=CB_color_cycle[2], fontsize=15)
    else:
        ax.plot(_wave[mask]/1e4, _throughput[mask], lw=2, color=CB_color_cycle[2], alpha=0.5, label=r"$\textbf{NIRES (order 3-7)}$")
        ax.plot(_wave[mask]/1e4, _throughput[mask] / _slit_loss[mask], lw=1, color='black', ls="-.", label="slit loss corrected")
        ax.plot(_wave[mask]/1e4, _throughput[mask] / _slit_loss[mask] * _slit_loss_psf[mask], lw=1, color='black', ls="dotted", label="seeing=0.6 arcsec")
        ax.text(_wave[mask][-1]/1e4, _throughput[mask][-1], f"{7-j}", color=CB_color_cycle[2], fontsize=15)

# ax.plot([], [], color="grey", lw=1, ls="dashed", label="Seeing = 0.9 arcsec for \nLRIS, MOSFIRE and \n NIRES order 7")
# ax.plot([], [], color="grey", lw=1, ls="-.", label="Slit loss corrected")

ax.set_ylim(-0.01, 0.3)
ax.set_xlabel(r'Wavelength ($\mu m$)', fontsize=20)
ax.set_ylabel('Throughput', fontsize=20)
ax.legend(fontsize=12, loc="lower right")
ax.tick_params(axis='both', which='major', labelsize=18)


wave_min, wave_max = ax.get_xlim()
wave_min, wave_max = wave_min*1e4, wave_max*1e4
z_min, z_max = wave_min / 1215.67 - 1, wave_max / 1215.67 - 1
ax2 = ax.twiny()
ax2.set_xlim(z_min, z_max)
ax2.set_xlabel(r"Ly$\alpha$ redshift", fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=18)

plt.show()