# Auto-generated PypeIt file using PypeIt version: 1.8.2.dev481+g93e899a45
# 2022-06-26

# User-defined execution parameters
[fluxcalib]
  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm

# Please add your SENSFUNC file name below before running pypeit_flux_calib

# Read in the flux
flux read
 Science_coadd/spec1d_r220305-r220305-J1100+0203.fits ../GD153_lris_sens.fits
 Science_coadd/spec1d_r220305-r220305-J1200+0112.fits
flux end

