# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-06-27

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_m220409-m220409-J1150.fits GD153_mosfire_sens.fits
 Science_coadd/spec1d_m220409-m220409-J1202+0129.fits
 Science_coadd/spec1d_m220409-m220409-J1315+1533.fits
 Science_coadd/spec1d_m220409-m220409-J1318+2932.fits
 Science_coadd/spec1d_m220409-m220409-J1332+0150.fits
 Science_coadd/spec1d_m220409-m220409-J1333+0919.fits
 Science_coadd/spec1d_m220409-m220409-J1511+0344.fits
 Science_coadd/spec1d_m220409-m220409-z7.fits
flux end
