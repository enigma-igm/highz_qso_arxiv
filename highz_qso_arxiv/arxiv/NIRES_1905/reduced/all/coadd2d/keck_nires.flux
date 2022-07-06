# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-07-06

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True
  
# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_NR.20190519.21141-NR.20190519.30380-J1007+2115.fits ../GD153_nires_sens.fits
 Science_coadd/spec1d_NR.20190519.33237-NR.20190519.33598-J1432+5746.fits
 Science_coadd/spec1d_NR.20190519.34377-NR.20190519.34740-J1458+6300.fits
 Science_coadd/spec1d_NR.20190519.35512-NR.20190519.35874-J1450+3302.fits
 Science_coadd/spec1d_NR.20190519.38953-NR.20190519.41150-J1535+1943.fits
 Science_coadd/spec1d_NR.20190519.45067-NR.20190519.48396-J1635+1758.fits
 Science_coadd/spec1d_NR.20190519.49989-NR.20190519.51855-J2212+2040.fits
flux end

