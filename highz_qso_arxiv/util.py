import imp
import numpy as np
from astropy.time import Time
from astropy.table import Table
from urllib.request import urlopen
from bs4 import BeautifulSoup

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

def get_skyprobe_extinction(t_start, t_end, format="unix"):
    """get extinction data from SkyProbe
       link: https://www.cfht.hawaii.edu/Instruments/Elixir/skyprobe/archive-new.html

        NOTE: there is a small, but variable, zero-point offset to the extinction measurement. 
              This varies based on sky conditions, instrument cleanliness (dust on the lens), 
              and telescope airmass. 
              For simplicity's sake, CFHT staff usually assume an average offset of +0.03.
              Thus we probably should -0.03 from the extinction we get.
    Args:
        t_start (int): UNIX time
        t_end (int):UNIX time
    """

    if format != "unix":
        if format != "mjd":
            raise TypeError("Time format unsupported")
        elif format == "mjd":
            t_start = mjd_to_unix(t_start)
            t_end = mjd_to_unix(t_end)

    url = f"https://www.cfht.hawaii.edu/ObsInfo/Weather/current/skyprobe_data.php?b={t_start}&e={t_end}"
    html = urlopen(url, timeout=10)
    soup = BeautifulSoup(html, 'html.parser')
    table_str = soup.get_text()

    # turn string into actual table (or list)
    table_str = table_str.replace("\n","")
    table_str = table_str.split("\r")
    table = np.array([r.split(",") for r in table_str])

    t_unix = table[:,0].astype("int") # UNIX time
    extinction = table[:,1].astype("float") # SkyProbe extinction measurement, in mag
    scatter = table[:,2].astype("float") # scatter of measurement, in mag

    # construct Table
    data = Table([t_unix, extinction, scatter], 
                 names=("Time", "Extinction", "Scatter"), 
                 units=("UNIX", "mag", "mag"))
    return data