import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from pypeit import io
from pypeit.sensfunc import IRSensFunc
from pypeit.core.flux_calib import Nlam_to_Flam

from IPython import embed

def plot_sens(sens=None, sens_file=None, title=None, display=True):
    if sens_file is not None:
        hdul = io.fits_open(sens_file)
        sens = IRSensFunc.from_hdu(hdul)

    wave = sens.wave
    zeropoint = sens.zeropoint
    throughput = sens.throughput
    sensitivity = (3631 * u.Jy*u.s*u.AA) * 10**(-zeropoint/2.5) * (c.c/(wave*u.AA)**2)
    sensitivity = sensitivity.to(u.erg/u.cm**2).value
    # sensitivity = Nlam_to_Flam(wave, zeropoint)

    fig, ax = plt.subplots(3, 1, figsize=(8, 8))
    ax[0].plot(wave, throughput)
    ax[0].set_xlabel('Wavelength [Angstrom]')
    ax[0].set_ylabel('Throughput')
    _, ymax = ax[0].get_ylim()
    ax[0].set_ylim(0, ymax)

    ax[1].plot(wave, zeropoint)
    ax[1].set_xlabel('Wavelength [Angstrom]')
    ax[1].set_ylabel('Zero Point [AB mag]')
    
    ax[2].plot(wave, sensitivity)
    ax[2].set_xlabel('Wavelength [Angstrom]')
    ax[2].set_ylabel('Sensitivity [erg/cm^2/photons]')

    fig.suptitle(title)
    fig.tight_layout()
    if display:
        plt.show()
    return fig, ax

if __name__ == "__main__":
    plot_sens(sens_file="../resource/GD153_mosfire_sens.fits", title="MOSFIRE")
    plot_sens(sens_file="../resource/GD153_lris_sens.fits", title="LRIS")