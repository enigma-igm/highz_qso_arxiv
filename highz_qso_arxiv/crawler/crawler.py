from highz_qso_arxiv.util import mjd_to_unix

from astropy.table import Table
from urllib.request import urlopen
from bs4 import BeautifulSoup
import numpy as np

from IPython import embed

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
                 names=("time", "extinction", "scatter"), 
                 units=("UNIX", "mag", "mag"))
    return data