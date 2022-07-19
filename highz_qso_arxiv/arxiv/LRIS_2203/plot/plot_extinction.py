import imp
from highz_qso_arxiv.util import mjd_to_unix
from highz_qso_arxiv.crawler import get_skyprobe_extinction
from pypeit.par.util import parse_pypeit_file
import matplotlib.pyplot as plt
import numpy as np

pypeit_file = '../LRIS_220305/reduced/all/keck_lris_red_mark4_Q.pypeit'
cfg_lines, data_files, frametype, usrdata, setups, _ = parse_pypeit_file(pypeit_file, runtime=True)
science_frame_mask = usrdata["frametype"] == "science"
obs_mjd = usrdata["mjd"][science_frame_mask]
obs_unix = [mjd_to_unix(mjd) for mjd in obs_mjd]

extinction = np.zeros(len(obs_mjd))
for i, unix in enumerate(obs_unix):
    try:
        extinction_data = get_skyprobe_extinction(unix-1800, unix+1800)
    except ValueError:
        raise ValueError("{}: No extinction data found for {}-{}".
                         format(usrdata["filename"][i], unix-1800, unix+1800))
    closest_time = min(extinction_data["time"], key=lambda x:abs(x-unix))
    dat = extinction_data[extinction_data["time"]==closest_time]
    extinction[i] = dat["extinction"].value[0]-0.03
