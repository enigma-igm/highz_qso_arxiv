import speclite.filters

ukirt_J_dat = ascii.read("../resource/UKIRT_UKIDSS.J.dat")
ukirt_J = speclite.filters.FilterResponse(
    wavelength=ukirt_J_dat['col1'] * u.AA,
    response=ukirt_J_dat["col2"], meta=dict(group_name='UKIRT', band_name='J')
)
ukirt = speclite.filters.load_filters('UKIRT-J')
