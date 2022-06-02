# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev86+g649552d6c
# 2022-04-19

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science/spec1d_m220409_0036-J0841+3814_OFF_MOSFIRE_20220409T064407.772.fits LDS749B_mosfire_sens.fits
 Science/spec1d_m220409_0037-J0841+3814_OFF_MOSFIRE_20220409T064710.362.fits
 Science/spec1d_m220409_0218-LDS749B_MOSFIRE_20220409T153212.032.fits
 Science/spec1d_m220409_0219-LDS749B_MOSFIRE_20220409T153444.012.fits
 Science_coadd/spec1d_m220409-m220409-J0841+3814.fits
flux end

