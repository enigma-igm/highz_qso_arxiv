# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-09-18

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm
  extrap_sens = True

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_MF.20201023.48683-MF.20201023.49610-J0339.fits ../GD71_mosfire_sens.fits
flux end

