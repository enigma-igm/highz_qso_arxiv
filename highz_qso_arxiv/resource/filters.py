import os
import speclite.filters
import astropy.units as u
from astropy.io import ascii

from highz_qso_arxiv.util import get_project_root

root = get_project_root()

ukirt_J_dat = ascii.read(os.path.join(root, 'resource/filter/UKIRT_UKIDSS.J.dat'))
ukirt_J = speclite.filters.FilterResponse(
    wavelength=ukirt_J_dat['col1'] * u.AA,
    response=ukirt_J_dat["col2"], meta=dict(group_name='UKIRT', band_name='J')
)
ukirt_J = speclite.filters.load_filters('UKIRT-J')
hsc_z = speclite.filters.load_filters('hsc2017-z')
decam_z = speclite.filters.load_filters('decam2014-z')
sdss_i = speclite.filters.load_filters('sdss2010-i')
wise_W1 = speclite.filters.load_filters('wise2010-W1')
