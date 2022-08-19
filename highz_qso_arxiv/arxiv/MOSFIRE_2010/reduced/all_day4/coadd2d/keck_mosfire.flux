# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-08-16

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True
  
# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_MF.20201025.17083-MF.20201025.18002-z7.fits ../GD71_mosfire_sens.fits
 Science_coadd/spec1d_MF.20201025.19017-MF.20201025.19939-z7.fits
 Science_coadd/spec1d_MF.20201025.22750-MF.20201025.22933-z7.fits
 Science_coadd/spec1d_MF.20201025.51880-MF.20201025.52804-z7.fits
flux end

