import imp
import numpy as np
from pathlib import Path
from astropy.time import Time

from IPython import embed

def inverse(array):
    return (array > 0.0)/(np.abs(array) + (array == 0.0))

def ivarsmooth(flux, ivar, window):
    """
    Boxcar smoothign of width window with ivar weights
    Args:
        flux:
        ivar:
        window:
    Returns:
    """

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
        smoothflux = np.roll(flux, 2*halfwindow + 1) # kill off NAN's

    return smoothflux, outivar

def luminosity_to_flux(luminosity, luminosity_distance):
    """
    convert luminosity to flux
    """
    return luminosity / (4 * np.pi * luminosity_distance**2)

def redshift_to_distance(redshift, cosmo):
    """
    Returns the distance in Mpc for a given redshift.
    """
    return cosmo.luminosity_distance(redshift).value
    
def mjd_to_unix(t_in_mjd):
    t = Time(t_in_mjd, format="mjd")
    return t.unix

def mjd_to_iso(t_in_mjd):
    t = Time(t_in_mjd, format="mjd")
    return t.iso

def unix_to_mjd(t_in_unix):
    t = Time(t_in_unix, format="unix")
    return t.mjd

def unix_to_iso(t_in_unix):
    t = Time(t_in_unix, format="unix")
    return t.iso

def get_project_root() -> Path:
    return Path(__file__).parent.parent
